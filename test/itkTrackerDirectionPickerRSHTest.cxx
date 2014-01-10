/**
 * @file  itkSymRealSphericalHarmonicRepTest.cxx
 * @brief Test itkSymRealSphericalHarmonicRep module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkTrackerDirectionPickerRSH.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbe.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>
#include <vector>

namespace{


// template <typename DirContainerType, typename DirectionType>
// int
// findClosestDirection(const DirContainerType allDirs, const DirectionType dir)
// {
//   double bestInner = 0.0;
//   double innerProd;
//   int bestIndex = -1;

//   for (unsigned int i = 0; i < allDirs.size(); ++i)
//   {
//     innerProd = allDirs[i][0]*dir[0] + allDirs[i][1]*dir[1] + allDirs[i][2]*dir[2];
//     if (bestInner < innerProd)
//     {
//       bestInner = innerProd;
//       bestIndex = i;
//     }  
//   }
//   // std::cerr << dir << std::endl;
//   // std::cerr << bestIndex << " : " << allDirs[bestIndex] << std::endl;
//   return bestIndex;
// }

// template <typename DirContainerType, typename DirectionType,typename HistType>
// DirectionType
// findMaxDirection(const DirContainerType allDirs, const HistType hist)
// {
//   double maxVal = 0.0;
//   int bestIndex = -1;

//   for (unsigned int i = 0; i < allDirs.size(); ++i)
//   {
//     if (maxVal < hist[i])
//     {
//       maxVal = hist[i];
//       bestIndex = i;
//     }  
//   }
//   // std::cerr << dir << std::endl;
//   // std::cerr << bestIndex << " : " << allDirs[bestIndex] << std::endl;
//   return allDirs[bestIndex];
// }


int testInitialize()
{

  const unsigned int OdfOrder = 4;
  typedef double                                                          PrecisionType;
  typedef itk::SymRealSphericalHarmonicRep<PrecisionType,OdfOrder>        PixelType;

  typedef itk::Image<PixelType,3>          ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                           InterpType;
  
  typedef itk::TrackerDirectionPickerRSH<InterpType>
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

  const unsigned int OdfOrder = 4;
  typedef double                                                          PrecisionType;
  typedef itk::SymRealSphericalHarmonicRep<PrecisionType,OdfOrder>        PixelType;

  typedef itk::Image<PixelType,3>           ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                            InterpType;

  typedef itk::TrackerDirectionPickerRSH<InterpType>
                                            PickerType;

  typedef ImageType::DirectionType          CosignMatrixType;

  CosignMatrixType LPSmat;
  LPSmat.Fill(0.0);
  LPSmat[0][0] = 1;
  LPSmat[1][1] = 1;
  LPSmat[2][2] = 1;

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

  double  c3[45] = {
               2.82094806e-01, -1.84499600e-03, -4.37358674e-03, -2.57216394e-01,
              -2.21597566e-03, -5.69619704e-03,  3.87160391e-01, -5.81498304e-03,
              -4.88071050e-03,  3.55852279e-03,  1.55444667e-01,  1.25033641e-03,
               2.96612363e-03,  2.67166062e-03, -6.83821831e-03, -1.60360534e-03,
              -3.04495008e-03, -1.28530160e-01,  2.44779466e-03,  5.92144812e-03,
              -4.25606035e-04, -6.48434013e-02,  5.08837402e-04, -1.10527303e-03,
              -2.69610246e-05,  2.88396259e-03, -2.25038873e-03, -7.32504530e-03,
               1.14771903e-01, -3.58927739e-03, -2.93563283e-03,  1.99614253e-04,
               2.79063080e-02,  2.29189231e-04, -2.89570610e-03, -1.19774416e-03,
               1.54047897e-02, -7.50655425e-04,  2.04480690e-04, -1.01385545e-03,
              -9.12652526e-04, -5.92372544e-06,  1.25703122e-03,  2.75207148e-03,
              -3.13866278e-03
            };

  PixelType odf1(c3);
  
  image->SetPixel(index,odf1);

  InterpType::Pointer interp  = InterpType::New();
  interp->SetInputImage(image);

  PickerType::Pointer pt = PickerType::New();
  pt->SetInterpolator(interp);

  pt->SetMethod(PickerType::STT);
  if (pt->GetMethod() != PickerType::STT)
  {
    std::cout << "[FAILED] setting picker method to STT failed" << std::endl;
    return EXIT_FAILURE;
  }
  pt->SetMethod(PickerType::STT);

  PickerType::DirectionType dir;
  PickerType::DirectionType expectedDir;
  PickerType::DirectionType lastDir;

 //TEST LAS ------------------------------------
  image->SetDirection(LPSmat);

  std::cout << "Testing Starting direction" << std::endl;
  dir = pt->PickStartingDirection(point);
  std::cout << dir << std::endl;

  std::cout << "Testing Pick Next direction" << std::endl;
  itk::TimeProbe * calc2Probe = new itk::TimeProbe(); 
  calc2Probe->Start();
  lastDir[0] = 1;  lastDir[1] = 0;  lastDir[2] = 0;
  std::cout << pt->PickNextDirection(point,lastDir) << std::endl;
  lastDir[0] = 0;  lastDir[1] = 1;  lastDir[2] = 0;
  std::cout << pt->PickNextDirection(point,lastDir) << std::endl;
  lastDir[0] = 0;  lastDir[1] = 0;  lastDir[2] = 1;
  std::cout << pt->PickNextDirection(point,lastDir) << std::endl;
  calc2Probe->Stop();
  std::cout << "PickNextDirection 3 elapsed time: " << calc2Probe->GetMean() << std::endl;

  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}

// int
// testPROB()
// {
//   std::cout << "Testing Probablistic direcction picking" << std::endl;

//   typedef itk::DiffusionTensor3D<double>    PixelType;
//   typedef itk::Image<PixelType,3>           ImageType;
//   typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
//                                             InterpType;

//   typedef itk::TrackerDirectionPickerDTI<InterpType>
//                                             PickerType;

//   typedef ImageType::DirectionType          CosignMatrixType;


//   typedef std::vector<PickerType::DirectionType> DirContainerType;
//   typedef vnl_vector<double>                HistType;

//   CosignMatrixType LPSmat;
//   LPSmat.Fill(0.0);
//   LPSmat[0][0] = 1;
//   LPSmat[1][1] = 1;
//   LPSmat[2][2] = 1;

//   CosignMatrixType PSLmat;
//   PSLmat.Fill(0.0);
//   PSLmat[0][2] = 1;
//   PSLmat[1][0] = 1;
//   PSLmat[2][1] = 1;

//   ImageType::RegionType region;
//   ImageType::IndexType start;
//   start[0] = 0;  start[1] = 0;  start[2] = 0;
 
//   ImageType::SizeType size;
//   size[0] = 1;  size[1] = 1;  size[2] = 1;
 
//   region.SetSize(size);
//   region.SetIndex(start);
 
//   ImageType::Pointer image = ImageType::New();
//   image->SetRegions(region);
//   image->Allocate();
  
//   ImageType::IndexType index;
//   index[0] = 0;  index[1] = 0;  index[2] = 0;

//   ImageType::PointType point;
//   image->TransformIndexToPhysicalPoint(index,point);

//   //this points straight up the z-axis....
//   PixelType dti;
//   dti[0] = 0.2;
//   dti[1] = 0;
//   dti[2] = 0;
//   dti[3] = 0.2;
//   dti[4] = 0;
//   dti[5] = 1.7;

//   image->SetPixel(index,dti);

//   InterpType::Pointer interp  = InterpType::New();
//   interp->SetInputImage(image);

//   PickerType::Pointer pt = PickerType::New();
//   pt->SetInterpolator(interp);
//   // pt->UseProbabilisticMethod();

//   PickerType::DirectionType dir;
//   PickerType::DirectionType expectedDir;
//   PickerType::DirectionType lastDir;
//   PickerType::DirectionType trueOdfPeak;


//   //We will test my making a histogram of selected directions

//   //1) make a bunch of gradient directions...
//   DirContainerType gradDirs;

//   unsigned int nDirs = 512;
//   double theta, phi, lastPhi=0.0, h;
//   for (unsigned int m = 0; m < nDirs; m++)
//   {
//     /*** Grad directions relation to theta phi
//     * x = sin(theta) * cos(phi)
//     * y = sin(theta) * sin(phi)
//     * z = cos(theta)
//     */
//     h = -1 + 2*( m )/ static_cast<double>( nDirs -1 );
//     theta = vcl_acos(h);
//     if (m == 0 or m == nDirs -1)
//     {
//       phi = 0.0;
//     }
//     else
//     {
//       phi = vcl_fmod( lastPhi + 3.6 / vcl_sqrt( nDirs * (1- vcl_pow(h,2) ) )
//                                         , ( 2 * vnl_math::pi ) );
//     }
//     dir[0] = sin(theta) * cos(phi); 
//     dir[1] = sin(theta) * sin(phi);
//     dir[2] = cos(theta);
//     gradDirs.push_back(dir);
//     lastPhi = phi;
//   }


//   HistType trueOdfHist(nDirs,0.0);
//   for (unsigned int i=0;i<nDirs;++i)
//   {
//     trueOdfHist[i] = EvaluateOdf(dti,gradDirs[i]);
//   }
//   trueOdfPeak = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,trueOdfHist);
//   std::cerr << trueOdfPeak << std::endl;

//   unsigned int N = nDirs * 100;

//   //TEST LPS ------------------------------------
//   std::cout << "Testing LPS direction set..." << std::endl;
//   image->SetDirection(LPSmat);

//   //generate LPS histogram
//   DirContainerType lpsDirs;
//   HistType lpsHist(nDirs,0.0);
//   for (unsigned int i=0;i<N;++i)
//   {
//     dir = pt->PickStartingDirection(point);
//     lpsHist[ findClosestDirection(gradDirs,dir) ] += 1.0/N;
//     lpsDirs.push_back(dir);
//   }
//   dir = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,lpsHist);
//   std::cerr << dir << std::endl;
//   //expectedDir[0] = 0;expectedDir[1] = 0;expectedDir[2] = 1;

//   expectedDir = LPSmat * trueOdfPeak;
//   if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.99)
//   {
//     std::cout << "  Failed StartingDirection -- expecting : "<<expectedDir << std::endl;
//     return EXIT_FAILURE;
//   }
 
//   //Test PSL
//   std::cout << "Testing PSL direction set..." << std::endl;
//   image->SetDirection(PSLmat);

//   DirContainerType pslDirs;
//   HistType pslHist(nDirs,0.0);
//   for (unsigned int i=0;i<N;++i)
//   {
//     dir = pt->PickStartingDirection(point);
//     pslHist[ findClosestDirection(gradDirs,dir) ] += 1.0/N;
//     pslDirs.push_back(dir);
//   }
//   dir = findMaxDirection<DirContainerType,PickerType::DirectionType,HistType>(gradDirs,pslHist);
//   std::cerr << dir << std::endl;
// //  expectedDir[0] = 1; expectedDir[1] = 0; expectedDir[2] = 0;
//   expectedDir = PSLmat * trueOdfPeak;
//   if ( vcl_abs( dir[0]*expectedDir[0]+dir[0]*expectedDir[1]+dir[2]*expectedDir[2] )   < 0.99)
//   {
//     std::cout << "  Failed StartingDirection -- expecting : "<<expectedDir << std::endl;
//     return EXIT_FAILURE;
//   }

//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   std::cout << "Testing Probabilistic succeeded" << std::endl << std::endl;

//   return EXIT_SUCCESS;
// }



} // end empty namespace

using namespace itk;
int itkTrackerDirectionPickerRSHTest( int argc, char * argv[] )
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
