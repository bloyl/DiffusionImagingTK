/**
 * @file  itkReplaceSpecialFunctions.cxx
 * @brief Some numerical helper functions.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkReplaceSpecialFunctions_cxx
#define __itkReplaceSpecialFunctions_cxx

#include <vnl/vnl_math.h>
#include <cmath>
#include <cstdlib>

#include <stdio.h>
#include <iostream>

#include <itkReplaceSpecialFunctions.h>
 
/**
 * Compute the binomial coefficient
 * 
 *    (  n )
 *    (  k )
 *
 * Using Comb6 method described in 
 * Binomial Coefficient Computation: Recursion or Iteration?
 *     Yannis Manolopoulos
 *     manolopo@delab.csd.auth.gr
 */
int
binomialCoeff( int n, int k )
{
  if (k==0 || n==k)
  {
    return 1;
  }

  if (n<0 or k<0 or k > n)
  {
    return 0;
  }

  int result = 1;
  int i;
  if (k < n-k)
  {
    for(i=n; i>= n-k+1; i--)
    {
      result *= i;
      result /= (n-i+1);
    }
  }
  else
  {
    for(i=n; i>= k+1; i--)
    {
      result *= i;
      result /= (n-i+1);
    }
  }
  return result;
}

/** Compute the associated legendre polynomial P_n^m(z)
 * 
 *  Use the following recursion relations from wolfram
 *  [1] (n-m) P_n^m(x) = x(2n-1) P_(n-1)^m(x) - (n+m-1) P_(n-2)^m(x). 
 * 
 *  and initializing it using
 *  [2] P_n^n(x)       = (-1)^n (2n-1)!! (1-x^2)^(n/2) 
 *  [3] P_(n+1)^n(x)   = x (2n+1)P_n^n(x). 
 * 
 *  note that (2n-1)!! is the Product of all odd integers <= (2n-1)
 *    so (2n-1)!! = 2n-1 * 2n-3 * 2n-5 * ... * 5 * 3 * 1
 * 
 * Extend the definition to negative m using the following identity
 *  [4] P_n^m(x)       = (-1)^m((n+m)!)/((n-m)!)  P_n^(-m)(x). 
 * 
 * This matched the definition in Wolfram mathematic 7 of LegendreP[n,m,z]
 **/

double
LegendreP( int n, int m, double x )
{
  //Handle the case of negative m
  if (m < 0)
  {
    //Determine (-1)^(m)
    int sign = 1;
    if ( vcl_abs(m) % 2 == 1) sign = -1;

    // Compute factorial(n+m) / factorial(n-m) 
    double f = 1.0;
    // m < 0 so (n-m) > (n+m)
    // (n+m)!/(n-m)! = Prod
    for (int i = (n-m); i>(n+m) ; i--)
    {
      f /= i;
    }
    
    return sign * f * LegendreP(n,-m,x);
  }
  
  if (m == 0 && n == 0)
  {
    return 1.0;
  }
  
  //Compute P_{n-2)_m
  //initialize with P_m^m = (-1)^m (2m-1)!! (1-x^2)^(m/2)
  double pm2 = 1;
  if ( m % 2 == 1 )
  {
    pm2 = -1;
  }

  //(-1)^m (2m-1)!!
  for ( int i=1 ; i<2*m ; i=i+2 )
  {
    pm2 *= i;
  }
  
  pm2 *= vcl_pow(vcl_sqrt( (1+x)*(1-x) ), m);

  if (m==n) return pm2;

  //Compute P_(m+1)^m(x)   = x (2m+1)P_m^m(x).
  double pm1 = x * (2 * m + 1) * pm2;

  if (n==m+1) return pm1;
  
  // Iterate (n-m) P_n^m(x) = x(2n-1) P_(n-1)^m(x) - (n+m-1) P_(n-2)^m(x).
  double pn=0;
  for (int nn = m+2; nn<=n; nn++)
  {
    pn = (x * (2 * nn -1) * pm1 - (nn+m-1) * pm2) / (nn - m);
    pm2 = pm1;
    pm1 = pn;
  }
  return pn;

}

/** Compute deriviatives of the associated legendre polynomial P_n^m(z)
 * 

  Forumlas from mathematica
In[16]:= CForm[D[LegendreP[n, m, z], {z, 1}] // FullSimplify]
Out[16]//CForm=
(-((1 + n)*z*LegendreP(n,m,z)) + (1 - m + n)*LegendreP(1 + n,m,z))/(-1 + Power(z,2))

In[15]:= CForm[D[LegendreP[n, m, z], {z, 2}] // FullSimplify]
Out[15]//CForm=
((1 + n + (1 + n)*(2 + n)*Power(z,2))*LegendreP(n,m,z) + (-1 + m - n)*((5 + 2*n)*z*LegendreP(1 + n,m,z) + (-2 + m - n)*LegendreP(2 + n,m,z)))/Power(-1 + Power(z,2),2)

 **/

void
LegendreP_derivatives( int n, int m, double z, double &first, double &second )
{
  double Pnmz      = LegendreP(n,m,z);
  double Pnplus1mz = LegendreP(n+1,m,z);
  double Pnplus2mz = LegendreP(n+2,m,z);

  first  = (-((1 + n)*z*Pnmz) + (1 - m + n)*Pnplus1mz)/(-1 + z*z);
  second = ((1 + n + (1 + n)*(2 + n)*z*z)*Pnmz + (-1 + m - n)*((5 + 2*n)*z*Pnplus1mz + (-2 + m - n)*Pnplus2mz))/(-1 + z*z)/(-1 + z*z);
}






#endif //__itkReplaceSpecialFunctions_cxx
