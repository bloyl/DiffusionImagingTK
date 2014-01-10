/**
 * @file  itkReplaceSpecialFunctions.h
 * @brief Some numerical helper functions.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkReplaceSpecialFunctions_h
#define __itkReplaceSpecialFunctions_h

/// Evaluates the Associated Legendre Polynomial with parameters (l,m) at x.
double LegendreP( int l, int m, double x );
///Compute the Binomial coefficient
int binomialCoeff( int n, int k );
///Compute the first and second derivatives of the legendre polynomials...
void LegendreP_derivatives( int n, int m, double z, double &first, double &second );

#endif //__itkReplaceSpecialFunctions_h
