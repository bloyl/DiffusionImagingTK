/**
* @file  itkDiffusionModelCalculator.hxx
* @brief This class is an abstract class that serves as the base for classes used
* to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkDiffusionModelCalculator_hxx
#define __itkDiffusionModelCalculator_hxx

#include <itkDiffusionModelCalculator.h>

#include <itkVector.h>

#include <vnl/algo/vnl_svd.h>
#include <vnl/vnl_rank.h>

namespace itk
{

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::DiffusionModelCalculator()
{
  m_NumGenerator = GeneratorType::New();
  m_NumGenerator->Initialize();

  m_GradientDirectionContainer = NULL;
  m_BMatrixInverse.set_identity();
  
  m_BeltramiLambda = 0.0; //By default use no Regularization.

  m_TensorReconMethod = DT_OLS;
  m_RSHReconMethod    = RSH_CSAODF;
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
void
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

///////////////////////////////////////////////////////////////////////////////////////////////

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
void
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::InitializeTensorFitting()
{

  if (!this->m_GradientDirectionContainer )
  {
    itkExceptionMacro( << "No gradient directions supplied. Need to supply at least 6" );
  }

  unsigned int numberOfGradientDirections = this->m_GradientDirectionContainer->Size();

  if( numberOfGradientDirections < 7 )
  {
    itkExceptionMacro( << "Not enough gradient directions supplied. Need to supply at least 7" );
  }

  m_BMatrix.set_size( numberOfGradientDirections, 7 );

  for (unsigned int m = 0; m < numberOfGradientDirections; m++)
  {
    m_BMatrix[m][0] = -     m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[0];
    m_BMatrix[m][1] = - 2 * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[1];
    m_BMatrix[m][2] = - 2 * m_GradientDirectionContainer->ElementAt(m)[0] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][3] = -     m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[1];
    m_BMatrix[m][4] = - 2 * m_GradientDirectionContainer->ElementAt(m)[1] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][5] = -     m_GradientDirectionContainer->ElementAt(m)[2] * m_GradientDirectionContainer->ElementAt(m)[2];
    m_BMatrix[m][6] = 1.0;
  }

  if (vnl_rank (m_BMatrix) < 7)
  {
    itkExceptionMacro("BMatrix is underdetermined -- more non colinear gradient direcitons are needed")
  }
  vnl_svd<double> solver(m_BMatrix);
  m_BMatrixInverse = solver.inverse();
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::DtType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ForceTensorToSPD(DtType dt) const
{         

  typename DtType::EigenValuesArrayType       lambda;
  typename DtType::EigenVectorsMatrixType     eigenVectors;

  typedef typename DtType::MatrixType         DTMatrixType;

  dt.ComputeEigenAnalysis(lambda, eigenVectors);
  bool tensorIsNotSpd = false;

  //Check all the eigenvalues
  for (unsigned int r=0; r<3; ++r)
  {
    if (lambda[r] <= 0)
    {
      tensorIsNotSpd = true;
      lambda[r] = -lambda[r];
    }
  }

  //if its notSPD then reconstitue the tensor..
  if (tensorIsNotSpd)
  {
    DTMatrixType diag;
    for (unsigned int r=0; r<3; ++r)
    {
      diag[r][r] = lambda[r];
    }

    DTMatrixType tmp
      = (static_cast<DTMatrixType>(eigenVectors.GetTranspose()) ) * (diag * eigenVectors);

    dt[0] = tmp[0][0];
    dt[1] = tmp[0][1];
    dt[2] = tmp[0][2];
    dt[3] = tmp[1][1];
    dt[4] = tmp[1][2];
    dt[5] = tmp[2][2];
  }
  return dt;
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::DtType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeTensor(DWIPixelType rawDWI, bool forceSPD) const
{
  DtType retVal;
  switch( m_TensorReconMethod )
  {
    case DT_OLS:
      retVal = ComputeTensorOLS(rawDWI,forceSPD);
      break;
    case DT_WLS:
      retVal = ComputeTensorWLS(rawDWI,forceSPD);
      break;
    case DT_WLSBOOT:
      retVal = ComputeTensorWLS_residualBoot(rawDWI,forceSPD);
      break;
    default:
      itkExceptionMacro( << "unrecognized Tensor Recon : " << m_TensorReconMethod << " : " << DT_OLS)
  }
  //Should never reach here.
  return retVal;
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::DtType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeTensorOLS(DWIPixelType rawDWI, bool forceSPD) const
{

  unsigned int nGrads = rawDWI.GetSize();

  //compute log of DWI signal
  vnl_vector<PrecisionType>  logDWI(nGrads, static_cast<PrecisionType>(0.0) );
  for (unsigned int i=0;i<nGrads;++i)
  {
    logDWI[i] = log(rawDWI[i]);
  }

  //Compute ordinary least squares data (mu_{OLS} between equations 9 and 10)
  vnl_vector<PrecisionType> olsDTVec =  m_BMatrixInverse * logDWI;

  // //Compute residual
  // vnl_vector<PrecisionType> olsFitData = m_BMatrix * olsDTVec;
  // for (int i=0;i<nGrads;++i)
  // {
  //   residual[i] = rawDWI[i] - exp(olsFitData[i]);
  // }

  DtType dt;
  dt(0,0) = olsDTVec[0];
  dt(0,1) = olsDTVec[1];
  dt(0,2) = olsDTVec[2];
  dt(1,1) = olsDTVec[3];
  dt(1,2) = olsDTVec[4];
  dt(2,2) = olsDTVec[5];

  if (forceSPD) 
    dt =  ForceTensorToSPD(dt);

  return dt;
}

//use weighted least squares (S Chung et al Neuroimage 33 (2006) p531-541)
template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::DtType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeTensorWLS(DWIPixelType rawDWI, bool forceSPD) const
{

  unsigned int nGrads = rawDWI.GetSize();

  //compute log of DWI signal
  vnl_vector<PrecisionType>  logDWI(nGrads, static_cast<PrecisionType>(0.0) );
  for (unsigned int i=0;i<nGrads;++i)
  {
    logDWI[i] = log(rawDWI[i]);
  }

  //Compute ordinary least squares data (mu_{OLS} between equations 9 and 10)
  vnl_vector<PrecisionType> olsFitData = m_BMatrix * m_BMatrixInverse * logDWI;

  //Compute the Weight matrix (defined between equations 9 and 10)
  CoefficientMatrixType W(nGrads,nGrads,0.0);
  for (unsigned int i=0;i<nGrads;++i)
  {
    PrecisionType tmp = exp(olsFitData[i]);
    W(i,i) = tmp * tmp;
  }

  //Compute part of the hatMatrix (around equation 11 and 12)
  // Hp = (X' W X)^(-1) X' W
  // H = X Hp
  vnl_svd<double> solver(m_BMatrix.transpose() * W * m_BMatrix);
  CoefficientMatrixType Hp = solver.inverse() * m_BMatrix.transpose() * W;
  CoefficientMatrixType H  = m_BMatrix * Hp;

  //Compute the weighted least square diffusion tensor (equation 9)
  vnl_vector<PrecisionType> wlsDt =  Hp * logDWI;

  DtType dt;
  dt(0,0) = wlsDt[0];
  dt(0,1) = wlsDt[1];
  dt(0,2) = wlsDt[2];
  dt(1,1) = wlsDt[3];
  dt(1,2) = wlsDt[4];
  dt(2,2) = wlsDt[5];

  if (forceSPD) 
    dt =  ForceTensorToSPD(dt);

  return dt;
}


template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::DtType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeTensorWLS_residualBoot(DWIPixelType rawDWI, bool forceSPD) const
{
  unsigned int nGrads = rawDWI.GetSize();

  //compute log of DWI signal
  vnl_vector<PrecisionType>  logDWI(nGrads, static_cast<PrecisionType>(0.0) );
  for (unsigned int i=0;i<nGrads;++i)
  {
    logDWI[i] = log(rawDWI[i]);
  }

  //Compute ordinary least squares data (mu_{OLS} between equations 9 and 10)
  vnl_vector<PrecisionType> olsFitData = m_BMatrix * m_BMatrixInverse * logDWI;

  //Compute the Weight matrix (defined between equations 9 and 10)
  CoefficientMatrixType W(nGrads,nGrads,0.0);
  for (unsigned int i=0;i<nGrads;++i)
  {
    PrecisionType tmp = exp(olsFitData[i]);
    W(i,i) = tmp * tmp;
  }

  //Compute part of the hatMatrix (around equation 11 and 12)
  // Hp = (X' W X)^(-1) X' W
  // H = X Hp
  vnl_svd<double> solver(m_BMatrix.transpose() * W * m_BMatrix);
  CoefficientMatrixType Hp = solver.inverse() * m_BMatrix.transpose() * W;
  CoefficientMatrixType H  = m_BMatrix * Hp;

  //Compute residuals for bootstrapping
  //scale each residual to have constant variance eq (11)
  vnl_vector<PrecisionType> residualErrors = logDWI - ( m_BMatrix * Hp * logDWI);

  PrecisionType meanRes = 0.0;
  for (unsigned int j=0;j<nGrads;++j)
  {
    residualErrors[j] *= sqrt(W(j,j));
    residualErrors[j] /= sqrt( 1 - H(j,j) );
    meanRes += residualErrors[j];
  }

  //Center by the mean
  meanRes /= nGrads;
  residualErrors -= meanRes;
  
  //Generate new dwiData...
  vnl_vector<PrecisionType>  logDWISample = olsFitData;
  for (unsigned int j=0;j<nGrads;++j)
  {
    int randIndex = m_NumGenerator->GetIntegerVariate(nGrads-1);
    logDWISample[j] += residualErrors[randIndex]/sqrt(W(j,j));
  }

  //REDO WLS fitting
  olsFitData = m_BMatrix * m_BMatrixInverse * logDWISample;

  //Don't think I need to zero W.
  //Compute the Weight matrix
  for (unsigned int i=0;i<nGrads;++i)
  {
    PrecisionType tmp = exp(olsFitData[i]);
    W(i,i) = tmp * tmp;
  }

  //Compute part of the hatMatrix
  // Hp = (X' W X)^(-1) X' W
  // H = X Hp
  vnl_svd<double> solver2(m_BMatrix.transpose() * W * m_BMatrix);
  Hp = solver2.inverse() * m_BMatrix.transpose() * W;
  
  vnl_vector<PrecisionType> wlsDtSamp =  Hp * logDWISample;

  DtType dt;
  dt(0,0) = wlsDtSamp[0];
  dt(0,1) = wlsDtSamp[1];
  dt(0,2) = wlsDtSamp[2];
  dt(1,1) = wlsDtSamp[3];
  dt(1,2) = wlsDtSamp[4];
  dt(2,2) = wlsDtSamp[5];

  if (forceSPD) 
    dt =  ForceTensorToSPD(dt);
  
  return dt;
}

///////////////////////////////////////////////////////////////////////////////////////////////

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
void
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::InitializeRSHFitting( TPrecisionType beltramiLambda = 0.006 )
{

  if (!this->m_GradientDirectionContainer )
  {
    itkExceptionMacro( << "No gradient directions supplied. Need to supply at least " << RSHNumberOfCoefficients );
  }

  unsigned int numberOfGradientDirections = this->m_GradientDirectionContainer->Size();

  if( numberOfGradientDirections < RSHNumberOfCoefficients)
  {
    itkExceptionMacro( << "Not enough gradient directions supplied. Need to supply at least " << RSHNumberOfCoefficients );
  }

  m_B0Indices.clear();
  m_DiffWeightedIndices.clear();

  //m_GradientDirectionContainer has dirs - (0,0,0) = b0 image as well as non unit length vectors
  // first find the B0 and DW Indices and normalize the gradient container
  typename GradientDirectionContainerType::Pointer DWIGradients = GradientDirectionContainerType::New();

  for(typename GradientDirectionContainerType::ConstIterator gdcit = m_GradientDirectionContainer->Begin();
      gdcit != m_GradientDirectionContainer->End(); ++gdcit)
  {
    GradientDirectionType gradDir = gdcit.Value();
    double gradNorm = gradDir.two_norm();
    if( gradNorm > 0.0 )
    {
      m_DiffWeightedIndices.push_back(gdcit.Index());
      DWIGradients->InsertElement(DWIGradients->Size(),gradDir / gradNorm);
    }
    else
    {
      m_B0Indices.push_back(gdcit.Index());      
    }
  }
  
  //Compute the forward basis for the DWIGradient directions
  m_RshBasis = RshType::ComputeRshBasis(DWIGradients);

  //compute the psuedo inverse once!
  if (beltramiLambda == 0)
  {
    vnl_svd<double> solver(m_RshBasis);
    m_RshBasisPseudoInverse = solver.inverse();
  }
  else
  {
    CoefficientMatrixType L2, basis2;
    typename RshType::LmVector lmVec;

    //Build lambda * L^2 matrix
    L2.set_size( RSHNumberOfCoefficients, RSHNumberOfCoefficients );
    L2.set_identity();

    for(unsigned int j=0;j<RSHNumberOfCoefficients;j++)
    {
      lmVec = RshType::GetLM(j+1);
      int l = lmVec(0);
      L2[j][j] = beltramiLambda * (l*l) * ( (l+1) * (l+1) );
    }
    basis2 = m_RshBasis.transpose() * m_RshBasis;
    vnl_svd<double> solver( basis2 + L2 );
    m_RshBasisPseudoInverse = solver.inverse() * m_RshBasis.transpose();
  }  
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::RshType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeRSH(DWIPixelType rawDWI) const
{
  RshType retVal;
  switch( m_RSHReconMethod )
  {
    case RSH_ADC:
      retVal = ComputeRSH_ADC(rawDWI);
      break;
    case RSH_ODF:
      retVal = ComputeRSH_ODF(rawDWI);
      break;
    case RSH_CSAODF:
      retVal = ComputeRSH_CSAODF(rawDWI);
      break;
    default:
      itkExceptionMacro( << "unrecognized RSH Recon")
  }
  return retVal;
}

template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::RshType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeRSH_ADC(DWIPixelType rawDWI) const
{
  //Compute s0
  typename NumericTraits<PrecisionType>::AccumulateType s0 = NumericTraits<PrecisionType>::Zero;
  for(unsigned int i = 0; i < m_B0Indices.size(); ++i)
  {
    s0 += rawDWI[m_B0Indices[i]];
  }
  s0 /= m_B0Indices.size();

  //Collect the log of the diffusion weighted single
  vnl_vector< PrecisionType > signal( m_DiffWeightedIndices.size() );

  for( unsigned int i=0; i < m_DiffWeightedIndices.size(); i++ )
  {
    if( rawDWI[m_DiffWeightedIndices[i]] == 0 )
    {
      signal[i] = 0;
    }
    else
    {
      signal[i] = log(static_cast<PrecisionType>( rawDWI[m_DiffWeightedIndices[i]] ) / static_cast<PrecisionType>(s0));
    }
    
    if ( vnl_math_isnan( rawDWI[m_DiffWeightedIndices[i]] )
          || vnl_math_isinf( rawDWI[m_DiffWeightedIndices[i]] ) )
    {
      itkExceptionMacro( "nans or Infs encountered\n" << rawDWI );
    }
  }

  vnl_vector< double > coeffs;
  coeffs  = m_RshBasisPseudoInverse * signal;
  
  if ( vnl_math_isnan(coeffs[1]) )
  {
    std::cerr << signal << std::endl;
    std::cerr << signal << std::endl;
    itkExceptionMacro( "nans or Infs encountered\n" << rawDWI );
  }

  RshType adc(0.0);
  for( unsigned int i=0; i < RSHNumberOfCoefficients; i++ )
    adc[i] = coeffs[i];  

  return adc;
}


template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::RshType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeRSH_ODF(DWIPixelType rawDWI) const
{
  //Compute s0
  typename NumericTraits<PrecisionType>::AccumulateType s0 = NumericTraits<PrecisionType>::Zero;
  for(unsigned int i = 0; i < m_B0Indices.size(); ++i)
  {
    s0 += rawDWI[m_B0Indices[i]];
  }
  s0 /= m_B0Indices.size();

  //Collect the diffusion weighted single
  vnl_vector< PrecisionType > signal( m_DiffWeightedIndices.size() );

  for( unsigned int i=0; i < m_DiffWeightedIndices.size(); i++ )
  {
    if( rawDWI[m_DiffWeightedIndices[i]] == 0 )
    {
      signal[i] = 0;
    }
    else
    {
      signal[i] = static_cast<PrecisionType>( rawDWI[m_DiffWeightedIndices[i]] ) / static_cast<PrecisionType>(s0);
    }
    if ( vnl_math_isnan( rawDWI[m_DiffWeightedIndices[i]] ) || vnl_math_isinf( rawDWI[m_DiffWeightedIndices[i]] ) )
    {
      itkExceptionMacro( "nans or Infs encountered\n" << rawDWI );
    }
  }

  vnl_vector< double > coeffs;
  coeffs  = m_RshBasisPseudoInverse * signal;
  
  //convert the rshcoeffs to odf coeffs
  for( unsigned int i=0; i<RSHNumberOfCoefficients; i++)
  {
    int l = (RshType::GetLM(i+1))[0];
    coeffs[i] = static_cast<double>( 2.0 * vnl_math::pi * LegendreP( l , 0, 0 ) * coeffs[i] );
  }

  //only normalize if 0th order is non zero...
  if ( coeffs[0] != 0)
  {
    double nFactor = 1.0 / ( 2 * vcl_sqrt( vnl_math::pi ) * coeffs[0]);
    for( unsigned int i=0; i< RSHNumberOfCoefficients; i++)
      coeffs[i] = nFactor * coeffs[i];
  }
  
  RshType odf(0.0);
  for( unsigned int i=0; i < RSHNumberOfCoefficients; i++ )
    odf[i] = static_cast<PrecisionType>(coeffs[i]);  

  return odf;
}


template < class TDWIPixelType, class TPrecisionType, unsigned int TOutputOrder >
typename DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>::RshType
DiffusionModelCalculator<TDWIPixelType, TPrecisionType, TOutputOrder>
::ComputeRSH_CSAODF( DWIPixelType rawDWI, TPrecisionType epsi1 = 0.15, TPrecisionType epsi2 = 0.15 ) const
{

  //Compute S0
  typename NumericTraits<PrecisionType>::AccumulateType s0 = NumericTraits<PrecisionType>::Zero;
  for(unsigned int i = 0; i < m_B0Indices.size(); ++i)
  {
    s0 += rawDWI[m_B0Indices[i]];
  }
  s0 /= m_B0Indices.size();

  //Compute E 
  vnl_vector< PrecisionType > signal( m_DiffWeightedIndices.size() );
  for( unsigned int i=0; i < m_DiffWeightedIndices.size(); i++ )
  {
    PrecisionType val = static_cast<PrecisionType>(rawDWI[m_DiffWeightedIndices[i]] / s0);
    PrecisionType regularizedVal;
    if ( val < 0 )
      regularizedVal = epsi1/2.0;
    else if ( val < epsi1 )
      regularizedVal = epsi1/2.0 + val * val / 2.0 / epsi1;
    else if ( val < (1-epsi2) )
      regularizedVal = val;
    else if ( val < 1 )
      regularizedVal = 1 - epsi2/2.0 - (1 - val)*(1 - val) / 2.0 / epsi2;
    else
      regularizedVal = 1 - epsi2/2.0;

    signal[i] = static_cast<PrecisionType>( log(-log( regularizedVal ) ) );
    
  }

  vnl_vector< PrecisionType > coeffs;
  coeffs  = m_RshBasisPseudoInverse * signal;

  // Apply 1/(16*pi^2) * LP to coeffs (ARound eq 3)
  for( unsigned int i=0; i<RSHNumberOfCoefficients; i++)
  {
    int l = (RshType::GetLM(i+1))[0];
    //  1 / (16 * pi^2) * -l * (l+1) * 2  pi * Legendre(l,0,0) = (-l * (l+1) * Pl(0,0) ) / 8 / pi
    coeffs[i] = static_cast<PrecisionType>( -l * (l+1) * LegendreP( l , 0, 0 ) * coeffs[i] / vnl_math::pi / 8.0 );
  }

  coeffs[0] += 1.0 / ( 2 * vcl_sqrt( vnl_math::pi ) );

  RshType odf(0.0);
  for( unsigned int i=0; i < RSHNumberOfCoefficients; i++ )
    odf[i] = static_cast<PrecisionType>(coeffs[i]);  

  return odf;
}



























}// end namespace itk

#endif
