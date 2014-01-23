
#include <itkSymRealSphericalHarmonicRep.h>
#include <itkDiffusionModelCalculator.h>
#include <itkPeakFindingCalculatorGrid.h>
#include <itkPeakFindingCalculator.h>

#include <itkTimeProbe.h>
#include <itkVariableLengthVector.h>
#include <vnl/vnl_random.h>

#include <itkTestingMacros.h>

#include "testingUtils.h"

using namespace DiffusionImagingTK_testing;

namespace{

template <typename DWIPixelType, typename GradientDirectionContainerType, typename DtType>
DWIPixelType
generateTestData(typename GradientDirectionContainerType::Pointer gradContainer,
                DtType dt, double s0, double noiseAmp)

{
  const unsigned int numberOfGradientImages = gradContainer->Size();
  typename GradientDirectionContainerType::Element dir;

  DWIPixelType dwi;
  dwi.SetSize(numberOfGradientImages);
  dwi.Fill(0);

  //Convert the dt to vnl mats
  vnl_matrix_fixed<double,3,3> T;
  T[0][0] = dt[0];
  T[0][1] = T[1][0] = dt[1];
  T[0][2] = T[2][0] = dt[2];
  T[1][1] = dt[3];
  T[1][2] = T[2][1] = dt[4];
  T[2][2] = dt[5];

  vnl_random randGen;
  for (unsigned int i=0;i<numberOfGradientImages;++i)
  {
    dir = gradContainer->ElementAt(i);
    dwi[i] = s0 * exp(-dot_product(dir, T * dir));
    if (noiseAmp > 0)
    {
      dwi[i] += abs(noiseAmp * randGen.normal());
    }
  }
  return dwi;
}

//TODO add test for SPD...itkDiffusionModelCalculatorTest

int diffusionTensorTest()
{
  typedef itk::VariableLengthVector<unsigned int>                 DWIPixelType;

  typedef itk::DiffusionModelCalculator<DWIPixelType, double>     DTCalculatorType;
  typedef DTCalculatorType::GradientDirectionContainerType        GradientDirectionContainerType;
  typedef DTCalculatorType::DtType                                DtType;
  typedef DTCalculatorType::RshType                               RshType;

  GradientDirectionContainerType::Pointer gradDirs = generateGradientDirections<GradientDirectionContainerType>(3);
  double bValue = 1000.0;
  double bValueSqrt = vcl_sqrt(bValue);

  std::cout << "Using " << gradDirs->Size() << " gradient Directions" << std::endl;
  
  //Gradient needs to be scaled by the sqrt of the bvalue...
  for (unsigned int i=0;i<gradDirs->Size();++i)
  {
    gradDirs->SetElement(i, bValueSqrt * gradDirs->ElementAt(i) );
  }

  DtType inpDt1;
  DtType inpDt2;
  DtType inpDt3;
  
  inpDt1[0] = 0.2E-3;
  inpDt1[1] = 0;
  inpDt1[2] = 0;
  inpDt1[3] = 0.5E-3;
  inpDt1[4] = 0;
  inpDt1[5] = 1.7E-3;
  inpDt2 = inpDt1;

  inpDt3[0] = -0.2E-3;
  inpDt3[1] = 0;
  inpDt3[2] = 0;
  inpDt3[3] = 0.5E-3;
  inpDt3[4] = 0;
  inpDt3[5] = 1.7E-3;

  //Test No Noise!
  DWIPixelType sig1 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt1, 5000., 0.);
  DWIPixelType sig2 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt2, 5000., 0.);
  DWIPixelType sig3 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt3, 5000., 0.);

  DWIPixelType res1;
  res1.SetSize(sig1.Size());

  DTCalculatorType::Pointer dtCalc = DTCalculatorType::New();
  EXERCISE_BASIC_OBJECT_METHODS(dtCalc, DTCalculatorType);

  TRY_EXPECT_EXCEPTION(dtCalc->InitializeTensorFitting());

  GradientDirectionContainerType::Pointer gradDirsSmall = GradientDirectionContainerType::New();
  gradDirsSmall->Reserve(4);
  for (unsigned int i=0;i<4;++i)
  {
    gradDirsSmall->InsertElement(i, gradDirs->ElementAt(i));
  }

  dtCalc->SetGradientDirectionContainer(gradDirsSmall);
  TRY_EXPECT_EXCEPTION(dtCalc->InitializeTensorFitting());

  dtCalc->SetGradientDirectionContainer(gradDirs);
  TEST_SET_GET(gradDirs, dtCalc->GetGradientDirectionContainer());

  dtCalc->InitializeTensorFitting();

  //Test ordinary least squares (OLS)
  DtType dt1    = dtCalc->ComputeTensorOLS(sig1, false);
  DtType dt2    = dtCalc->ComputeTensorOLS(sig2, false);
  DtType dt3    = dtCalc->ComputeTensorOLS(sig3, false);
  DtType dt3spd = dtCalc->ComputeTensorOLS(sig3, true);

  std::cout << "dt ols fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares (WLS)
  dt1    = dtCalc->ComputeTensorWLS(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS(sig3, true);

  std::cout << "dt WLS fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares with residual boot straping(WLS)
  dt1    = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS_residualBoot(sig3, true);

  std::cout << "dt WLS residual boot strap fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  // Test with noise
  sig1 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt1, 5000., 500.);
  sig2 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt2, 5000., 500.);
  sig3 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt3, 5000., 500.);

  //Test ordinary least squares (OLS)
  dt1    = dtCalc->ComputeTensorOLS(sig1, false);
  dt2    = dtCalc->ComputeTensorOLS(sig2, false);
  dt3    = dtCalc->ComputeTensorOLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorOLS(sig3, true);

  std::cout << "dt ols fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares (WLS)
  dt1    = dtCalc->ComputeTensorWLS(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS(sig3, true);

  std::cout << "dt WLS fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares with residual boot straping(WLS)
  dt1    = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS_residualBoot(sig3, true);

  std::cout << "dt WLS residual boot strap fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;


  // dtCalc->InitializeRSHFitting();

  // RshType adc1  = dtCalc->ComputeRSH_ADC(sig1);
  // RshType odf1  = dtCalc->ComputeRSH_ODF(sig1);
  // RshType csaodf = dtCalc->ComputeRSH_CSAODF(sig1);
  // RshType csaodf_noReg = dtCalc->ComputeRSH_CSAODF(sig1,0,0);

  // std::cout << "adc1         " << adc1 << std::endl;
  // std::cout << "odf1         " << odf1 << std::endl;
  // std::cout << "csaodf       " << csaodf << std::endl;
  // std::cout << "csaodf_noReg " << csaodf_noReg << std::endl;

  // //Test Some Noise!
  // generateTestData<DWIPixelType, GradientDirectionContainerType> (gradDirs, bValue, inpDt1, inpDt2, inpDt3,5000,500, sig1, sig2, sig3);

  // dt1 = dtCalc->ComputeTensorOLS(sig1, false);
  // dt2 = dtCalc->ComputeTensorOLS(sig2, false);
  // dt3 = dtCalc->ComputeTensorOLS(sig3, false);

  // wlsDt1 = dtCalc->ComputeTensorWLS(sig1, false);
  // wlsDt2 = dtCalc->ComputeTensorWLS(sig2, false);
  // wlsDt3 = dtCalc->ComputeTensorWLS(sig3, false);

  // resBootDt1 = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  // resBootDt2 = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  // resBootDt3 = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);

  // std::cout << "sig1 " << std::endl;
  // std::cout << "ols         " << dt1 << std::endl;
  // std::cout << "wls         " << wlsDt1 << std::endl;
  // std::cout << "resBootDt   " << resBootDt1 << std::endl;

  // std::cout << "sig2 " << std::endl;
  // std::cout << "ols         " << dt2 << std::endl;
  // std::cout << "wls         " << wlsDt2 << std::endl;
  // std::cout << "resBootDt   " << resBootDt2 << std::endl;

  // std::cout << "sig3 " << std::endl;
  // std::cout << "ols         " << dt3 << std::endl;
  // std::cout << "wls         " << wlsDt3 << std::endl;
  // std::cout << "resBootDt   " << resBootDt3 << std::endl;

  // //Do some timing tests...
  // itk::TimeProbe * olsProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * wlsProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * resBootProbe = new itk::TimeProbe(); 
  // olsProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorOLS(sig1, false);
  // olsProbe->Stop();

  // wlsProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorWLS(sig1, false);
  // wlsProbe->Stop();

  // resBootProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  // resBootProbe->Stop();

  // std::cout << "200 ols elapsed time            : " << olsProbe->GetMean() << std::endl;
  // std::cout << "200 wls elapsed time            : " << wlsProbe->GetMean() << std::endl;
  // std::cout << "200 res Bootstrap elapsed time  : " << resBootProbe->GetMean() << std::endl;

  // std::cout << "*************************************************************************************" << std::endl
  //           << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
  //           << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}

} // end empty namespace

using namespace itk;
int itkDiffusionModelCalculatorTest( int , char ** )
{

  if ( diffusionTensorTest() == EXIT_FAILURE )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
