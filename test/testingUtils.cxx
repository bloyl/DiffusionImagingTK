
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