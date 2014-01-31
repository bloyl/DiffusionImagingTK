/*=========================================================================
 *
 *  Copyright Section of Biomedical Image Analysis
 *            University of Pennsylvania
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include <itkTrackerDirectionPickerDWI.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkImage.h>
#include <itkVectorImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>
#include <vector>

namespace{

// int testInitialize()
// {
//   typedef itk::VectorImage<double, 3>      GradientImageType;
//   typedef itk::NearestNeighborInterpolateImageFunction<GradientImageType,double>
//                                            InterpType;
  
//   typedef itk::TrackerDirectionPickerDWI<InterpType>
//                                            PickerType;
//   typedef PickerType::Superclass           PickerBaseType;
  
//   PickerType::Pointer pt = PickerType::New();
//   std::cerr << pt << std::endl;

//   PickerBaseType::Pointer spt = dynamic_cast<PickerBaseType *>(pt.GetPointer());
 
//   std::cerr << spt << std::endl;

//   return EXIT_FAILURE;
// }

///////////////////////////////////////////////////////////////////////////////////////////

template <typename GradientImageType, typename GradientDirectionContainerType>
void
generateGradientImage(typename GradientImageType::Pointer gradientImage,
   typename GradientDirectionContainerType::Pointer gradContainer,
   double bValue)
{
  const unsigned int numberOfGradientImages = 67; // The bbl set!!!
  // const unsigned int numberOfGradientImages = 8; // The bbl set!!!

  //Make 3 Signals
  unsigned int sig1[numberOfGradientImages] =
  {
    5566, 6309, 4317, 1797, 3370, 2046, 6236, 3713, 3346, 5616, 4226,
       5315, 2358, 3962, 4287, 3713, 7129, 6478, 2757, 4000, 6020, 5192,
       4384, 5929, 4569, 1789, 4600, 3977, 3023, 3764, 4483, 2498, 1404,
       1183, 2655, 4127, 5258, 5463, 5191, 3770, 7130, 2400, 3632, 2150,
       3442, 2484, 2959, 1248, 3911, 5410, 2275, 6601, 3598, 2230, 3907,
       3104, 2762, 5018, 3489, 3050, 2960, 1682, 2509, 3082, 1467, 4418,
       3488
  };

  unsigned int sig2[numberOfGradientImages] =
  {
    5403, 1174, 3430, 3990, 4854, 4320, 2467, 1588, 5105, 3344, 3080,
       4076, 3926, 4780, 3039, 2969, 2717, 3212, 4396, 4304, 3727, 3743,
       2501, 2789, 4345, 6994, 3380, 3325, 4287, 4168, 2992, 3968, 4684,
       4665, 4407, 2080, 7469, 1158, 2673, 3143, 2127, 2285, 3333, 6163,
       2797, 3415, 2964, 4651, 2793, 2214, 4169, 4370, 3008, 3795, 4357,
       4823, 4724, 5562, 3972, 4354, 3650, 4986, 3455, 6089, 4426, 2348,
       4467
  };

  unsigned int sig3[numberOfGradientImages] =
  {
    5566, 6309, 4317, 1797, 3370, 2046, 6236, 3713, 3346, 5616, 4226,
       5315, 2358, 3962, 4287, 3713, 7129, 6478, 2757, 4000, 6020, 5192,
       4384, 5929, 4569, 1789, 4600, 3977, 3023, 3764, 4483, 2498, 1404,
       1183, 2655, 4127, 5258, 5463, 5191, 3770, 7130, 2400, 3632, 2150,
       3442, 2484, 2959, 1248, 3911, 5410, 2275, 6601, 3598, 2230, 3907,
       3104, 2762, 5018, 3489, 3050, 2960, 1682, 2509, 3082, 1467, 4418,
       3488
  };

  typedef typename GradientImageType::Pointer              GradientImagePointer;
  typedef typename GradientImageType::RegionType           GradientRegionType;
  typedef typename GradientRegionType::IndexType           GradientIndexType;
  typedef typename GradientRegionType::SizeType            GradientSizeType;
  typedef typename GradientImageType::PixelType            GradPixelType;

  typedef typename GradientImageType::PointType            PointType;


  gradientImage->SetVectorLength(numberOfGradientImages);
  
  GradientSizeType  sizeGradientImage  = {{ 1, 1, 3 }};
  GradientIndexType indexGradientImage = {{ 0, 0, 0 }};
  GradientRegionType     regionGradientImage;
  regionGradientImage.SetSize(  sizeGradientImage );
  regionGradientImage.SetIndex( indexGradientImage);
  gradientImage->SetRegions( regionGradientImage );
  gradientImage->Allocate();

  PointType point;

  itk::ImageRegionIteratorWithIndex< GradientImageType > git(
        gradientImage, regionGradientImage );

  git.GoToBegin();

  /// lets fill in the signal where we want to...
  ///Odf0 in {0,0,0}
  GradPixelType pix(numberOfGradientImages);
  pix.Fill(0);
  for (unsigned int i=0;i<numberOfGradientImages;++i)
    pix[i] = sig1[i];
  indexGradientImage[0] = 0;
  indexGradientImage[1] = 0;
  indexGradientImage[2] = 0;
  gradientImage->SetPixel(indexGradientImage,pix);
  gradientImage->TransformIndexToPhysicalPoint(indexGradientImage,point);
  std::cerr << "Point 1 : " << point << std::endl;


  ///Odf1 in {0,0,1}
  indexGradientImage[0] = 0;
  indexGradientImage[1] = 0;
  indexGradientImage[2] = 1;
  pix.Fill(0);
  for (unsigned int i=0;i<numberOfGradientImages;++i)
    pix[i] = sig2[i];
  gradientImage->SetPixel(indexGradientImage,pix);
  gradientImage->TransformIndexToPhysicalPoint(indexGradientImage,point);
  std::cerr << "Point 2 : " <<  point << std::endl;

  ///Odf2 in {0,0,2}
  indexGradientImage[0] = 0;
  indexGradientImage[1] = 0;
  indexGradientImage[2] = 2;
  pix.Fill(0);
  for (unsigned int i=0;i<numberOfGradientImages;++i)
    pix[i] = sig3[i];
  gradientImage->SetPixel(indexGradientImage,pix);
  gradientImage->TransformIndexToPhysicalPoint(indexGradientImage,point);
  std::cerr << "Point 3 : " <<  point << std::endl;

  //Set up Gradient Contatiner...
  typedef typename GradientDirectionContainerType::Element       GradientDirectionType;
  
  GradientDirectionType dir;
  double  gradientDirections[numberOfGradientImages][3] =
  {
    {0.000000,0.000000,0.000000},
    {1.000000,0.000000,0.000000},
    {0.000000,1.000000,0.000000},
    {-0.026007,0.649170,0.760199},
    {0.591136,-0.766176,0.252058},
    {-0.236071,-0.524158,0.818247},
    {-0.893021,-0.259006,0.368008},
    {0.796184,0.129030,0.591137},
    {0.233964,0.929855,0.283956},
    {0.935686,0.139953,0.323891},
    {0.505827,-0.844710,-0.174940},
    {0.346220,-0.847539,-0.402256},
    {0.456968,-0.630956,-0.626956},
    {-0.486997,-0.388997,0.781995},
    {-0.617845,0.672831,0.406898},
    {-0.576984,-0.104997,-0.809978},
    {-0.826695,-0.520808,0.212921},
    {0.893712,-0.039987,-0.446856},
    {0.290101,-0.541189,-0.789276},
    {0.115951,-0.962591,-0.244896},
    {-0.800182,0.403092,-0.444101},
    {0.513981,0.839970,0.173994},
    {-0.788548,0.152912,-0.595659},
    {0.949280,-0.233069,0.211062},
    {0.232964,0.782880,0.576911},
    {-0.020999,-0.187990,-0.981946},
    {0.216932,-0.955701,0.198938},
    {0.774003,-0.604002,0.190001},
    {-0.160928,0.355840,0.920587},
    {-0.147035,0.731173,-0.666158},
    {0.888141,0.417066,0.193031},
    {-0.561971,0.231988,-0.793959},
    {-0.380809,0.142928,0.913541},
    {-0.306000,-0.199000,-0.931001},
    {-0.332086,-0.130034,0.934243},
    {-0.963226,-0.265062,0.044010},
    {0,0,0},
    {-0.959501,0.205107,0.193101},
    {0.452965,-0.888932,0.067995},
    {-0.773133,0.628108,0.088015},
    {0.709082,0.408047,0.575066},
    {-0.692769,0.023992,0.720760},
    {0.681659,0.528735,-0.505747},
    {-0.141995,-0.724976,0.673978},
    {-0.740168,0.388088,0.549125},
    {-0.103006,0.822044,0.560030},
    {0.584037,-0.596038,0.551035},
    {-0.088008,-0.335031,0.938088},
    {-0.552263,-0.792377,0.259123},
    {0.838158,-0.458086,-0.296056},
    {0.362995,-0.560993,0.743990},
    {-0.184062,0.392133,-0.901306},
    {-0.720938,-0.692941,0.008999},
    {0.433101,0.682159,-0.589137},
    {0.502114,0.690157,0.521119},
    {-0.170944,-0.508833,-0.843722},
    {0.462968,0.422971,0.778946},
    {0,0,0},
    {0.385030,-0.809064,0.444035},
    {-0.713102,-0.247035,0.656094},
    {0.259923,0.884737,-0.386885},
    {0.001000,0.077002,-0.997030},
    {0.037002,-0.902057,0.430027},
    {0.570320,-0.303170,-0.763428},
    {-0.282105,0.145054,-0.948354},
    {0.721098,0.608082,0.332045},
    {0.266985,0.959945,-0.084995}
  };
  
  double BValueSqrt = vcl_sqrt(bValue);
  for (unsigned int g = 0; g<numberOfGradientImages;++g)
  {
    dir[0] = BValueSqrt * gradientDirections[g][0];
    dir[1] = BValueSqrt * gradientDirections[g][1];
    dir[2] = BValueSqrt * gradientDirections[g][2];
    gradContainer->InsertElement(g,dir);
  }
}


int test1()
{
  typedef itk::VectorImage<double, 3>       GradientImageType;
  typedef GradientImageType::Pointer        GradientImagePointer;
  typedef GradientImageType::PointType      PointType;

  typedef itk::NearestNeighborInterpolateImageFunction<GradientImageType,double>
                                            InterpType;
  
  typedef itk::TrackerDirectionPickerDWI<InterpType>
                                            PickerType;
  typedef PickerType::Superclass            PickerBaseType;
  
  typedef PickerType::GradientDirectionContainerType  GradientDirectionContainerType;


  //Make the gradient Image...
  GradientImagePointer gradIm = GradientImageType::New();
  GradientDirectionContainerType::Pointer gradDirs = GradientDirectionContainerType::New();
  generateGradientImage<GradientImageType,GradientDirectionContainerType>(gradIm,gradDirs,1000.0);

  InterpType::Pointer interp  = InterpType::New();
  interp->SetInputImage(gradIm);

  PickerType::Pointer pt = PickerType::New();
  pt->SetInterpolator(interp);
  pt->SetGradientDirectionContainer(gradDirs);

  pt->Initialize();
  PointType point(0.0);
  PickerType::DirectionType dir;

  point[0] = 0;  point[1] = 0;  point[2] = 0;
  dir = pt->PickStartingDirection(point);
  std::cerr << point << " : " << dir << std::endl;
  
  point[0] = 0;  point[1] = 0;  point[2] = 1;
  dir = pt->PickStartingDirection(point);
  std::cerr << point << " : " << dir << std::endl;

  point[0] = 0;  point[1] = 0;  point[2] = 2;
  dir = pt->PickStartingDirection(point);
  std::cerr << point << " : " << dir << std::endl;

  pt->UseCSAODFStreamLineTracking();

  point[0] = 0;  point[1] = 0;  point[2] = 0;
  dir = pt->PickStartingDirection(point);
  std::cerr << point << " : " << dir << std::endl;
  // dir = pt->PickNextDirection(point,dir);
  // std::cerr << point << " :N " << dir << std::endl;

  // point[0] = 0;  point[1] = 0;  point[2] = 1;
  // dir = pt->PickStartingDirection(point);
  // std::cerr << point << " : " << dir << std::endl;
  // dir = pt->PickNextDirection(point,dir);
  // std::cerr << point << " :N " << dir << std::endl;

  // point[0] = 0;  point[1] = 0;  point[2] = 2;
  // dir = pt->PickStartingDirection(point);
  // std::cerr << point << " : " << dir << std::endl;
  // dir = pt->PickNextDirection(point,dir);
  // std::cerr << point << " :N " << dir << std::endl;

  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;
            
  return EXIT_FAILURE;
}


} // end empty namespace

using namespace itk;
int itkTrackerDirectionPickerDWITest( int , char ** )
{
  // if ( testInitialize() )
  //   return EXIT_FAILURE;

  if ( test1() )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
