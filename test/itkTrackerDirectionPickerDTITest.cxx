/**
 * @file  itkSymRealSphericalHarmonicRepTest.cxx
 * @brief Test itkSymRealSphericalHarmonicRep module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkTrackerDirectionPickerDTI.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>
#include <vector>


namespace{


template <typename DirContainerType, typename DirectionType>
int
findClosestDirection(const DirContainerType allDirs, const DirectionType dir)
{
  double bestInner = 0.0;
  double innerProd;
  int bestIndex = -1;

  for (unsigned int i = 0; i < allDirs.size(); ++i)
  {
    innerProd = allDirs[i][0]*dir[0] + allDirs[i][1]*dir[1] + allDirs[i][2]*dir[2];
    if (bestInner < innerProd)
    {
      bestInner = innerProd;
      bestIndex = i;
    }  
  }
  // std::cerr << dir << std::endl;
  // std::cerr << bestIndex << " : " << allDirs[bestIndex] << std::endl;
  return bestIndex;
}

template <typename DirContainerType, typename DirectionType,typename HistType>
DirectionType
findMaxDirection(const DirContainerType allDirs, const HistType hist)
{
  double maxVal = 0.0;
  int bestIndex = -1;

  for (unsigned int i = 0; i < allDirs.size(); ++i)
  {
    if (maxVal < hist[i])
    {
      maxVal = hist[i];
      bestIndex = i;
    }  
  }
  // std::cerr << dir << std::endl;
  // std::cerr << bestIndex << " : " << allDirs[bestIndex] << std::endl;
  return allDirs[bestIndex];
}


int testInitialize()
{
  typedef itk::DiffusionTensor3D<double>   PixelType;
  typedef itk::Image<PixelType,3>          ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                           InterpType;
  
  typedef itk::TrackerDirectionPickerDTI<InterpType>
                                           PickerType;
  typedef PickerType::Superclass           PickerBaseType;
  
  PickerType::Pointer pt = PickerType::New();
  std::cerr << pt << std::endl;

  PickerBaseType::Pointer spt = dynamic_cast<PickerBaseType *>(pt.GetPointer());
 
  std::cerr << spt << std::endl;

  return EXIT_FAILURE;
}

int
testDeterminisitic()
{
  std::cout << "Testing Dterminisitic direction picking" << std::endl;
  typedef itk::DiffusionTensor3D<double>    PixelType;
  typedef itk::Image<PixelType,3>           ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                            InterpType;

  typedef itk::TrackerDirectionPickerDTI<InterpType>
                                            PickerType;

  typedef ImageType::DirectionType          CosignMatrixType;

  CosignMatrixType LPSmat;
  LPSmat.Fill(0.0);
  LPSmat[0][0] = 1;
  LPSmat[1][1] = 1;
  LPSmat[2][2] = 1;

  CosignMatrixType PSLmat;
  PSLmat.Fill(0.0);
  PSLmat[0][2] = 1;
  PSLmat[1][0] = 1;
  PSLmat[2][1] = 1;

  std::cerr << PSLmat << std::endl;

  ImageType::RegionType region;
  ImageType::IndexType start;
  start[0] = 0;  start[1] = 0;  start[2] = 0;
 
  ImageType::SizeType size;
  size[0] = 1;  size[1] = 1;  size[2] = 1;
 
  region.SetSize(size);
  region.SetIndex(start);
 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  
  ImageType::IndexType index;
  index[0] = 0;  index[1] = 0;  index[2] = 0;

  ImageType::PointType point;
  image->TransformIndexToPhysicalPoint(index,point);

  //this points straight up the z-axis....
  PixelType dti;
  dti[0] = 0.2;
  dti[1] = 0;
  dti[2] = 0;
  dti[3] = 0.5;
  dti[4] = 0;
  dti[5] = 1.7;

  image->SetPixel(index,dti);

  InterpType::Pointer interp  = InterpType::New();
  interp->SetInputImage(image);

  PickerType::Pointer pt = PickerType::New();
  pt->SetInterpolator(interp);

  PickerType::DirectionType dir;
  PickerType::DirectionType expectedDir;
  PickerType::DirectionType lastDir;

  //TEST LAS ------------------------------------
  image->SetDirection(LPSmat);

  std::cout << "Testing Starting direction" << std::endl;
  dir = pt->PickStartingDirection(point);
  expectedDir[0] = 0;expectedDir[1] = 0;expectedDir[2] = 1;
  if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  ---Failed StartingDirection(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  ---Passed StartingDirection(LPS) " << std::endl;

  //Test PSL
  image->SetDirection(PSLmat);

  dir = pt->PickStartingDirection(point);
  expectedDir[0] = 1; expectedDir[1] = 0; expectedDir[2] = 0;
  if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  ---Failed StartingDirection(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  ---Passed StartingDirection(PSL) " << std::endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << std::endl << "Testing NextDirection Picking" << std::endl;

  std::cout << "  ---Testing Streamline tracking" << std::endl;
  pt->UseStreamLineTracking();
  image->SetDirection(LPSmat);
  lastDir[0] = 0;  lastDir[1] = 0;  lastDir[2] = -1;
  expectedDir[0] = 0; expectedDir[1] = 0; expectedDir[2] = -1;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Streamline PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Streamline PickNext(LPS)" << std::endl;


  image->SetDirection(PSLmat);
  lastDir[0] = -1;  lastDir[1] = 0;  lastDir[2] = 0;
  expectedDir[0] = -1;  expectedDir[1] = 0; expectedDir[2] = 0;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Streamline PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Streamline PickNext(PSL)" << std::endl;



  std::cout << "  ---Testing Tensor Deflection tracking" << std::endl;
  pt->UseTensorDeflectionTracking();
  image->SetDirection(LPSmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = 0.7920;
  expectedDir[0] = 0.0214807; expectedDir[1] = -0.214807; expectedDir[2] = 0.97642;

  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensor Deflection PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensor Deflection PickNext(LPS)" << std::endl;


  image->SetDirection(PSLmat);
  lastDir[0] = 0.7920;  lastDir[1] = 0.1481;  lastDir[2] = -0.5924;
  expectedDir[0] = 0.97642; expectedDir[1] = 0.0214807; expectedDir[2] = -0.214807;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensor Deflection PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensor Deflection PickNext(PSL)" << std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "  ---Testing Tensorline tracking" << std::endl;
  pt->UseTensorLineTracking();

  pt->SetTensorLineF(1);
  pt->SetTensorLineG(0);

  std::cout << "  ------Testing f = " << pt->GetTensorLineF() << std::endl;
  std::cout << "  ------Testing g = " << pt->GetTensorLineG() << std::endl;

  image->SetDirection(LPSmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = -0.7920;
  expectedDir[0] = 0; expectedDir[1] = 0; expectedDir[2] = -1;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(LPS)" << std::endl;

  image->SetDirection(PSLmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = 0.7920;
  expectedDir[0] = 1;  expectedDir[1] = 0; expectedDir[2] = 0;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(PSL)" << std::endl;

  pt->SetTensorLineF(0);
  pt->SetTensorLineG(1);

  std::cout << "  ------Testing f = " << pt->GetTensorLineF() << std::endl;
  std::cout << "  ------Testing g = " << pt->GetTensorLineG() << std::endl;

  image->SetDirection(LPSmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = 0.7920;
  expectedDir[0] = 0.0214807; expectedDir[1] = -0.214807; expectedDir[2] = 0.97642;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(LPS)" << std::endl;

  image->SetDirection(PSLmat);
  lastDir[0] = 0.7920;  lastDir[1] = 0.1481;  lastDir[2] = -0.5924;
  expectedDir[0] = 0.97642; expectedDir[1] = 0.0214807; expectedDir[2] = -0.214807;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(PSL)" << std::endl;


  pt->SetTensorLineF(0);
  pt->SetTensorLineG(0.7);

  std::cout << "  ------Testing f = " << pt->GetTensorLineF() << std::endl;
  std::cout << "  ------Testing g = " << pt->GetTensorLineG() << std::endl;

  image->SetDirection(LPSmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = 0.7920;
  expectedDir[0] = 0.0607057; expectedDir[1] = -0.334921; expectedDir[2] = 0.940288;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(LPS)" << std::endl;

  image->SetDirection(PSLmat);
  lastDir[0] = 0.7920;  lastDir[1] = 0.1481;  lastDir[2] = -0.5924;
  expectedDir[0] = 0.940288; expectedDir[1] = 0.0607057; expectedDir[2] = -0.334921;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(PSL)" << std::endl;


  pt->SetTensorLineF(0.3);
  pt->SetTensorLineG(0.7);

  std::cout << "  ------Testing f = " << pt->GetTensorLineF() << std::endl;
  std::cout << "  ------Testing g = " << pt->GetTensorLineG() << std::endl;

  image->SetDirection(LPSmat);
  lastDir[0] = 0.1481;  lastDir[1] = -0.5924;  lastDir[2] = 0.7920;
  expectedDir[0] = 0.043037; expectedDir[1] = -0.237441; expectedDir[2] = 0.970448;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(LPS)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(LPS)" << std::endl;

  image->SetDirection(PSLmat);
  lastDir[0] = 0.7920;  lastDir[1] = 0.1481;  lastDir[2] = -0.5924;
  expectedDir[0] = 0.970448; expectedDir[1] = 0.043037; expectedDir[2] = -0.237441;
  dir = pt->PickNextDirection(point,lastDir);
  if ( ( dir[0]*expectedDir[0]+dir[1]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.9999)
  {
    std::cout << "  -----Failed Tensorline PickNext(PSL)" << std::endl;
    std::cout << "  ---expecting : "<<expectedDir << std::endl;
    std::cout << "  ---got       : "<<dir << std::endl;
    return EXIT_FAILURE;
  }
  else
    std::cout << "  -----Passed Tensorline PickNext(PSL)" << std::endl;

  std::cout << "Testing Determinisitic methods succeeded" << std::endl << std::endl;

  return EXIT_SUCCESS;
}

template <typename PixelType, typename DirectionType >
double
EvaluateOdf(const PixelType dti, const DirectionType dir)
{
  vnl_matrix_fixed<double,3,3> D;
  D(0,0)          = dti[0];
  D(0,1) = D(1,0) = dti[1];
  D(0,2) = D(2,0) = dti[2];
  D(1,1)          = dti[3];
  D(1,2) = D(2,1) = dti[4];
  D(2,2)          = dti[5];

  vnl_matrix_fixed<double,3,3> Dinv = vnl_inverse(D);
  const vnl_vector<double> g = dir.GetVnlVector();

  // val = 1 / (4 pi) * 1 / (sqrt(det(D)) * 1 / (g' D^(-1) * g)^(3/2)
  double val;
  val = 1. / 4. / vnl_math::pi;
  val /= vcl_sqrt( vnl_det(D) );
  val /= vcl_pow( inner_product(g,Dinv*g), 1.5);

  return val;
}

int
testPROB()
{
  std::cout << "Testing Probablistic direcction picking" << std::endl;

  typedef itk::DiffusionTensor3D<double>    PixelType;
  typedef itk::Image<PixelType,3>           ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                            InterpType;

  typedef itk::TrackerDirectionPickerDTI<InterpType>
                                            PickerType;

  typedef ImageType::DirectionType          CosignMatrixType;


  typedef std::vector<PickerType::DirectionType> DirContainerType;
  typedef vnl_vector<double>                HistType;

  CosignMatrixType LPSmat;
  LPSmat.Fill(0.0);
  LPSmat[0][0] = 1;
  LPSmat[1][1] = 1;
  LPSmat[2][2] = 1;

  CosignMatrixType PSLmat;
  PSLmat.Fill(0.0);
  PSLmat[0][2] = 1;
  PSLmat[1][0] = 1;
  PSLmat[2][1] = 1;

  ImageType::RegionType region;
  ImageType::IndexType start;
  start[0] = 0;  start[1] = 0;  start[2] = 0;
 
  ImageType::SizeType size;
  size[0] = 1;  size[1] = 1;  size[2] = 1;
 
  region.SetSize(size);
  region.SetIndex(start);
 
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  
  ImageType::IndexType index;
  index[0] = 0;  index[1] = 0;  index[2] = 0;

  ImageType::PointType point;
  image->TransformIndexToPhysicalPoint(index,point);

  //this points straight up the z-axis....
  PixelType dti;
  dti[0] = 0.2;
  dti[1] = 0;
  dti[2] = 0;
  dti[3] = 0.2;
  dti[4] = 0;
  dti[5] = 1.7;

  image->SetPixel(index,dti);

  InterpType::Pointer interp  = InterpType::New();
  interp->SetInputImage(image);

  PickerType::Pointer pt = PickerType::New();
  pt->SetInterpolator(interp);
  // pt->UseProbabilisticMethod();

  PickerType::DirectionType dir;
  PickerType::DirectionType expectedDir;
  PickerType::DirectionType lastDir;
  PickerType::DirectionType trueOdfPeak;


  //We will test my making a histogram of selected directions

  //1) make a bunch of gradient directions...
  DirContainerType gradDirs;

  unsigned int nDirs = 512;
  double theta, phi, lastPhi=0.0, h;
  for (unsigned int m = 0; m < nDirs; m++)
  {
    /*** Grad directions relation to theta phi
    * x = sin(theta) * cos(phi)
    * y = sin(theta) * sin(phi)
    * z = cos(theta)
    */
    h = -1 + 2*( m )/ static_cast<double>( nDirs -1 );
    theta = vcl_acos(h);
    if (m == 0 or m == nDirs -1)
    {
      phi = 0.0;
    }
    else
    {
      phi = vcl_fmod( lastPhi + 3.6 / vcl_sqrt( nDirs * (1- vcl_pow(h,2) ) )
                                        , ( 2 * vnl_math::pi ) );
    }
    dir[0] = sin(theta) * cos(phi); 
    dir[1] = sin(theta) * sin(phi);
    dir[2] = cos(theta);
    gradDirs.push_back(dir);
    lastPhi = phi;
  }


  HistType trueOdfHist(nDirs,0.0);
  for (unsigned int i=0;i<nDirs;++i)
  {
    trueOdfHist[i] = EvaluateOdf(dti,gradDirs[i]);
  }
  trueOdfPeak = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,trueOdfHist);
  std::cerr << trueOdfPeak << std::endl;

  unsigned int N = nDirs * 100;

  //TEST LPS ------------------------------------
  std::cout << "Testing LPS direction set..." << std::endl;
  image->SetDirection(LPSmat);

  //generate LPS histogram
  DirContainerType lpsDirs;
  HistType lpsHist(nDirs,0.0);
  for (unsigned int i=0;i<N;++i)
  {
    dir = pt->PickStartingDirection(point);
    lpsHist[ findClosestDirection(gradDirs,dir) ] += 1.0/N;
    lpsDirs.push_back(dir);
  }
  dir = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,lpsHist);
  std::cerr << dir << std::endl;
  //expectedDir[0] = 0;expectedDir[1] = 0;expectedDir[2] = 1;

  expectedDir = LPSmat * trueOdfPeak;
  if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.99)
  {
    std::cout << "  Failed StartingDirection -- expecting : "<<expectedDir << std::endl;
    return EXIT_FAILURE;
  }
 
  //Test PSL
  std::cout << "Testing PSL direction set..." << std::endl;
  image->SetDirection(PSLmat);

  DirContainerType pslDirs;
  HistType pslHist(nDirs,0.0);
  for (unsigned int i=0;i<N;++i)
  {
    dir = pt->PickStartingDirection(point);
    pslHist[ findClosestDirection(gradDirs,dir) ] += 1.0/N;
    pslDirs.push_back(dir);
  }
  dir = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,pslHist);
  std::cerr << dir << std::endl;
//  expectedDir[0] = 1; expectedDir[1] = 0; expectedDir[2] = 0;
  expectedDir = PSLmat * trueOdfPeak;
  if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.99)
  {
    std::cout << "  Failed StartingDirection -- expecting : "<<expectedDir << std::endl;
    return EXIT_FAILURE;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  std::cout << "Testing Probabilistic succeeded" << std::endl << std::endl;

  return EXIT_SUCCESS;
}



} // end empty namespace

using namespace itk;
int itkTrackerDirectionPickerDTITest( int, char ** )
{

  // if ( testInitialize() )
  //   return EXIT_FAILURE;
  if ( testDeterminisitic() == EXIT_FAILURE )
    return EXIT_FAILURE;
  
  // std::cout << "Testing Probabilistic Method:" << std::endl; 
  // if ( testPROB() == EXIT_FAILURE )
  //   return EXIT_FAILURE;


  return EXIT_SUCCESS;
}
