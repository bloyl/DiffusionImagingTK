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
