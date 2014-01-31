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


#ifndef __testingUtils_h
#define __testingUtils_h

namespace DiffusionImagingTK_testing
{

template <typename GradientDirectionContainerType>
void loadGradDirs( typename GradientDirectionContainerType::Pointer gradDirs, std::string bvecFile, std::string bvalFile );

template < class DwiImageType4D, class DwiImageType >
void convertDWIimage(typename DwiImageType4D::Pointer inp4d, typename DwiImageType::Pointer outputIm );

template <class GradientDirectionContainerType>
typename GradientDirectionContainerType::Pointer generateGradientDirections(int resolution);

bool areEqual(double x, double y, double percision);

}//End Namespace

//empty
int testingUtils(int, char **);

#endif

#include "testingUtils.hxx"
