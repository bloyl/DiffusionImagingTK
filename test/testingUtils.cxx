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


#ifndef __testingUtils_cxx
#define __testingUtils_cxx

#include <limits.h>
#include <vcl_cmath.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>      // std::setprecision

#include "testingUtils.h"

int testingUtils(int, char **)
{
  return EXIT_FAILURE;
}

namespace DiffusionImagingTK_testing
{

bool areEqual(double x, double y, double percision)
{
  //Compare binanry significant and eponent
  // x = xSig * 2^xExp
  double xSig,ySig;
  int xExp,yExp;

  xSig = frexp(x , &xExp);
  ySig = frexp(y , &yExp);

  if (xExp != yExp)
  {
    std::cerr << "areEqual exponents differ : " << xExp << " : " << yExp << std::endl;
    return false;
  }
  if ( vcl_abs(xSig - ySig) > percision )
  {
    return false;
  }
  return true;
}

}

#endif