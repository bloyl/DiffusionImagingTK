/**
 * @file  itkSymRealSphericalHarmonicRepTest.cxx
 * @brief Test itkSymRealSphericalHarmonicRep module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkTestingMacros.h>

#include <itkTrackerDirectionPickerImageBase.h>

#include <itkDiffusionTensor3D.h>
#include <itkNearestNeighborInterpolateImageFunction.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <limits.h>
#include <cmath>
#include <stdio.h>
#include <vector>


using namespace itk;
int itkTrackerDirectionPickerImageBaseTest( int, char ** )
{

  typedef itk::DiffusionTensor3D<double>   PixelType;
  typedef itk::Image<PixelType,3>          ImageType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double>
                                           InterpType;
  
  typedef itk::TrackerDirectionPickerImageBase<InterpType>
                                           PickerType;
  
  PickerType::Pointer picker = PickerType::New();
  EXERCISE_BASIC_OBJECT_METHODS( picker, PickerType );

  if (picker->IsInitialized())
  {
    std::cout << "picker is not intialized but returns that it is";
    return EXIT_FAILURE;
  }


  // InterpType::Pointer interp  = InterpType::New();
  // interp->SetInputImage(gradIm);

  // std::cout << "Passed" << std::endl;

  return EXIT_SUCCESS;
}
