/**
 * @file  itkPeakFindingCalculator.hxx
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPeakFindingCalculator_hxx
#define __itkPeakFindingCalculator_hxx


#include <itkPeakFindingCalculator.h>


#define  NALPLUT_ELERES  512    // number of elevations covered by the NALPs LUT
#define  NUM_SAMPLEDIRS  100    // number of guess directions --- samples on a unit sphere
#define  DIR_ARRAY_SIZE  200    // = NUM_SAMPLEDIRS * 2 --- two values per direction
#define  MAX_DIR_CHANGE  0.20
#define  MIN_DIR_CHANGE  1e-4
#define  MAX_DOTPRODUCT  0.99   // two fibers with the angle's absolute dot-product greater 
                                // than this value are considered co-linear, either duplicates
                                // or paired


// a set of guess directions, each by (phi, theta),
// sampled on a HEMI-sphere
static const double  SAMPLE_DIRS[ DIR_ARRAY_SIZE ] =   
       {  
          0,         0,            0,          1.05629, 
          2.22439,   0.545639,     0.432648,   1.38198,
         -2.66275,   0.459002,    -0.00929288, 1.65375,	
          2.65675,   0.644856,    -1.74803,    0.707892,	
         -2.27218,   1.266110,    -1.47540,    1.411200,	
          0.212308,  1.23103,      0.386965,   0.690593,	
         -0.947022,  1.4182,      -1.18591,    0.947845,	
         -0.859285,  1.66178,     -1.63763,    1.66476,	
          1.6568,    1.24269,      1.58155,    0.757079,	
          3.06922,   0.516202,     2.69357,    0.903883,	
          1.93779,   1.24967,      1.16823,    0.521205,	
          1.96264,   0.754403,     0.0308025,  0.795651,	
          1.20719,   0.781226,     2.42966,    1.05056,	
         -1.0762,    0.260592,    -1.21142,    1.40092,	
          2.11327,   0.988479,     2.97089,    1.03569,	
          1.24368,   1.50987,     -1.34743,    0.715331,	
         -0.445539,  1.36629,     -2.80037,    1.15395,	
         -2.19948,   1.65047,     -1.34293,    1.17168,	
          1.50272,   1.01189,     -1.50652,    0.948898,	
         -1.06013,   1.18294,      0.572587,   1.15978,	
          2.71186,   1.16407,     -0.238037,   1.52724,	
         -0.505823,  0.441938,    -1.11808,    1.64428,	
         -2.47615,   0.700818,    -2.90468,    1.39363,	
         -0.240991,  1.19776,     -1.58921,    0.459988,	
          0.291254,  0.977912,    -1.82455,    0.964427,	
         -0.27807,   0.940309,    -2.6669,     1.50328,	
          3.11708,   0.258538,     1.20234,    1.0431,	
          2.33938,   0.798317,    -1.62392,    1.18891,	
         -0.931793,  0.772514,    -0.513095,   1.11261,	
          3.00591,   0.775179,    -1.04543,    0.521355,	
         -2.116,     1.04884,      2.11867,    0.288387,	
         -0.0316375, 1.36282,      0.0580415,  0.259231,	
         -0.707266,  1.30266,     -0.625203,   1.55302,	
          0.805254,  0.68481,     -2.88577,    0.696834,	
          1.37843,   1.2482,      -2.1131,     0.789698,	
         -1.3776,    1.65648,      2.48143,    1.3084,	
          0.588921,  0.451236,    -2.12533,    0.529549,	
         -1.90086,   1.21922,      0.195789,   1.49202,	
          0.923059,  0.930546,    -2.06203,    1.42824,	
          1.80762,   1.0161,      -3.08051,    1.19604,	
         -2.40581,   1.49027,     -0.322049,   0.682883,	
          0.598741,  0.901629,    -0.595579,   0.860811,	
         -1.7387,    1.42551,      1.12527,    0.261047,	
         -0.806462,  1.05662,     -2.11988,    0.268098,	
         -0.191187,  1.84387,      1.10734,    1.28761,	
          0.0397623, 0.519667,    -2.68789,    0.913252,	
         -0.421495,  1.71465,      0.696674,   1.39439,	
         -2.41527,   1.03986,      1.70489,    0.504068,	
         -0.924502,  1.91203,      0.854966,   1.18342,	
         -2.54944,   1.26955,     -3.01216,    0.942817
       };


namespace itk
{


template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::PeakFindingCalculator()
{
  beReady = 0;
  toCheck = 1;
  LUTused = 1;
  nalpLUT = NULL;
  shCoefs.Fill( 0.0 );
  sortedPeaks.clear();

  rangeOfInterest[0] = 1.0e10;
  rangeOfInterest[1] =-1.0e10;

  GenerateNalpLUT();
}


template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::~PeakFindingCalculator()
{ 
  sortedPeaks.clear();

  if ( nalpLUT )
  {
    delete [] nalpLUT;
    nalpLUT = NULL;
  }
}


// Create an LUT for accelerated computation of NALP values
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
void PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::GenerateNalpLUT()
{
  if ( nalpLUT )   return;

  nalpLUT  =  new  TPrecisionType  [ Dimension * NALPLUT_ELERES ];
  
  TPrecisionType   eleVal = 0.0;
  TPrecisionType * elePtr = nalpLUT;
  TPrecisionType   eleInc = vnl_math::pi / ( NALPLUT_ELERES - 1 );

  for ( int i =  0; i < NALPLUT_ELERES; i ++, eleVal += eleInc, elePtr += Dimension ) 
  for ( int l =  0; l <= MaxOrder; l += 2 )
  for ( int m = -l; m <= l;        m ++   )
  {
    elePtr[ ShIndex( l, m ) ] = NALP( l, m, eleVal );
  }
  
  elePtr = NULL;
}


// set up the context for interpolating NALP values
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
void PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::SetInterpolator( TPrecisionType eleVal, Interp * interp )
{
  int               eleIdx;
  TPrecisionType    eleInc = vnl_math::pi / ( NALPLUT_ELERES - 1 );
  interp->weight1 = eleVal / eleInc;
  eleIdx = int( interp->weight1 );
  
  if ( eleIdx < 0 )
  {
    eleIdx = 0;
    interp->weight1 = 0.0;
    interp->weight0 = 1.0;
  }
  else
  if ( eleIdx >= NALPLUT_ELERES - 1 )
  {
    eleIdx = NALPLUT_ELERES - 1;
    interp->weight1 = 0.0;
    interp->weight0 = 1.0;
  }
  else
  {
    interp->weight1-= eleIdx;
    interp->weight0 = 1.0 - interp->weight1;
  }

  interp->lutPtr0 = nalpLUT + eleIdx * Dimension;
  interp->lutPtr1 = interp->lutPtr0  + Dimension;
}


// Set a range (max followed by min) of interest (in terms of the OD value) 
// within which to access the peaks.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
void PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::SetRangeOfInterest( OutputComponentType maxMin[2] )
{  
  rangeOfInterest[0] = maxMin[0];  
  rangeOfInterest[1] = maxMin[1];
  if ( maxMin[0] < maxMin[1] )
  {
    rangeOfInterest[0] = maxMin[1];
    rangeOfInterest[1] = maxMin[0];
  }
}


// The Normalization factor for spherical harmonic (l, m). xxx
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
TPrecisionType PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::K( int l, int m )
{
  PrecisionType  f = 1.0; // if m = 0

  if ( m > 0 )
  {
    for( int i = l - m + 1; i < l + m + 1; i ++ ) f /= i;
  }
  else
  {
    for( int i = l + m + 1; i < l - m + 1; i ++ ) f *= i;
  }

  return static_cast < PrecisionType >
         (    vcl_sqrt(  ( 2 * l + 1 ) * f / ( 4 * vnl_math::pi )  )    );
}


// This function computes the Normalized Associated Legendre Polynomial (NALP) of
// spherical harmonic function with order l and phase factor m for an elevation.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
TPrecisionType PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::NALP( int l, int m, PrecisionType elevation ) 
{
  if( m == 0 )
  {
    return static_cast < PrecisionType >
           (    PixelType::K( l, 0 ) 
              * LegendreP(  l,  m,  vcl_cos( elevation )  )    );
  }
  else // for m >0 or m < 0
  {
    // NOTE: The wiki-based spherical harmonic equation is a normalized version
    //       AND is complex-based. Thus the real and imaginary parts each have
    //       magnitude (1 / sqrt2). Since we will use the REAL part ONLY, we then
    //       need to scale up the wiki-based equation by sqrt2 to make magnitude 1.
    //
    return static_cast < PrecisionType >
           (    vnl_math::sqrt2 * PixelType::K( l, m )
              * LegendreP(  l,  m,  vcl_cos( elevation )  )    );
  }
}


// Given a total of ( MaxOrder + 1 ) * ( MaxOrder + 2 ) / 2 spherical harmonic
// coefficients as the representation an OD (either FOD or ODF), this function 
// computes the OD as well as 1st-order and 2nd-order derivatives along a specific
// direction (elevation, azimuth)
// 
// ODValue:     the OD value
// dSH_del:     the first-order  derivative of the OD relative to elevation
// dSH_daz:     the first-order  derivative of the OD relative to azimuth
// d2SH_del2:   the second-order derivative of the OD relative to elevation
// d2SH_deldaz: the second-order derivative of the OD relative to elevation & azimuth
// d2SH_daz2:   the second-order derivative of the OD relative to azimuth
//
// by linearly summing the values of each set of functions (represented via spherical
// harmonics, i.e., the OD set and the OD-derivative sets), with each function value
// weighted by the corresponding spherical harmonic coefficient. The OD and derivatives
// are then exploited for finding peaks within a voxel by means of the Newton-Raphson
// Gradient Descent method.
//
// As opposed to Tournier's representation, cos is assigned to negative phases
// and sin to positive phases in class itkSymRealSphericalHarmonicRep, which
// results in a different set of equations for calculating derivatives and OD
// values.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
void PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::GetODValueAndDerivatives( PrecisionType   elevation,   PrecisionType   azimuth,
                            PrecisionType * ODValue,     PrecisionType * dSH_del,
                            PrecisionType * dSH_daz,     PrecisionType * d2SH_del2,
                            PrecisionType * d2SH_deldaz, PrecisionType * d2SH_daz2 )
{
  int             l, m;
  int             nNALPs = ( MaxOrder + 1 ) * ( MaxOrder + 2 ) / 2; // ALL phase factors of EVEN orders
  PrecisionType   tmpVal = 0.0;
  PrecisionType   sinEle = static_cast < PrecisionType > (  vcl_sin( elevation )  );
  int             atPole = ( sinEle < 1e-4 ) ? 1 : 0;
  PrecisionType * pNALPs = new PrecisionType [ nNALPs ];
  
                       
  // init the values to be returned
  *ODValue = *dSH_del = *dSH_daz = *d2SH_del2 = *d2SH_deldaz = *d2SH_daz2 = 0.0;

  
  // compute a complete set of NALP values to be intensively reused
  if ( LUTused ) // accelerated computation through an LUT
  {
    Interp     interp;
    SetInterpolator( elevation, &interp );

    for ( l =  0; l <= MaxOrder; l += 2 )
    for ( m = -l; m <= l;        m ++   )
    {
      pNALPs[ ShIndex( l, m ) ] = NALP( &interp, l, m );
    }

    interp.lutPtr0 = interp.lutPtr1 = NULL;
  }
  else           // brute-force computation
  {
    for ( l =  0; l <= MaxOrder; l += 2 )
    for ( m = -l; m <= l;        m ++   )
    {
      pNALPs[ ShIndex( l, m ) ] = NALP( l, m, elevation );
    }
  }
  

  // 'pNALPs' stores the NALP values of the spherical harmonic functions at a given orientation
  // (elevation, azimuth). These values, each weighted by the sphrical harmonic coefficient
  // (provided by 'shCoefs'), are linearly combined together to obtain the overall OD value.
  //
  // function value = SUM{   ( spherical harmonic #i's coeficient     )
  //                       * ( spherical harmonic #i's function value )  }


  // Q( l, m, theta ) = sqrt(  ( 2 * l + 1 ) * ( l - m )! / 4PI / ( l + m )!  )
  //                  * P( l, m, theta )  
  //                  = NALP( l, m, theta ) --- simply denoted below as Q( l, m )
  //
  // Y( l, m, theta, phi ) = Q( l, m, theta ) * exp( i * m * phi )
  //
 

  // -------------------------------------------------------------------------------
  // the FIRST-order derivative of OD ---- dSH( l, m, theta, phi ) / dTheta
  // 
  // A1 = sqrt(  ( l + m ) * ( l - m + 1 )  )                B1 = Q( l, m - 1 )
  //
  // A2 = sqrt(  ( l - m ) * ( l + m + 1 )  )                B2 = Q( l, m + 1 )
  //
  // dQ( l, m ) / dTheta = ( A1 * B1 - A2 * B2 ) / 2
  //
  // dSH( l, m, theta, phi ) / dTheta =
  //		
  //    (1) m < 0     cos( m * phi ) * dQ( l, m ) / dTheta
  //
  //    (2) m > 0     sin( m * phi ) * dQ( l, m ) / dTheta
  //
  // NOTE: the cos-sin scheme here is just opposed to Tournier's


  // -------------------------------------------------------------------------------
  // the SECOND-order derivative of OD ---- dSH2( l, m, theta, phi ) / dTheta2
  //
  // A1 = sqrt(  ( l + m ) * ( l - m + 1 ) * ( l + m - 1 ) * ( l - m + 2 )  )
  //
  // B1 = Q( l, m - 2 )
  //
  // A2 = sqrt(  ( l - m ) * ( l + m + 1 ) * ( l - m - 1 ) * ( l + m + 2 )  )
  // B2 = Q( l, m + 2 )
  //
  // A3 = [ ( l + m ) * ( l - m + 1 ) + ( l - m ) * ( l + m + 1 ) ]
  // B3 = Q( l, m )
  //
  // dSH2( l, m, theta, phi ) / dTheta2 = ( A1 * B1 + A2 * B2 - A3 * B3 ) / 4


  // init three values with phase factor m = 0 --- a special treatment
  //
  for ( l = 0; l <= MaxOrder; l += 2 )
  { 
    // Any derivative relative to azimuth involves an m item and hence is ZERO for phase factor
    // 0. Thus in this initialization loop, we do not need to consider dSH_daz, dSH_deldaz, and
    // d2SH_daz2. Instead we only need to consider 3 values: ODValue, dSH_del, and d2SH_del2.

    // the OD value
    ( *ODValue ) += static_cast < PrecisionType > (  shCoefs[ ShIndex( l, 0 ) ]  )
                  * pNALPs[ ShIndex( l, 0 ) ];

    // dSH_del and d2SH_del2 --- derivatives relative to theta --- for non-zero orders ONLY
    //
    if ( l > 0 ) 
    {
      // the FIRST-order derivative of the OD
      //
      // A1 = A2 = sqrt(  l * ( l + 1 )  )
      //
      // B1 = Q( l, -1 )
      // B2 = Q( l,  1 )
      //
      // total = ( A1 * B1 - A2 * B2 ) / 2 = A1 * ( B1 - B2 ) / 2
      //
      ( *dSH_del ) -= static_cast < PrecisionType > (  shCoefs[ ShIndex( l, 0 ) ]  )     // xxx -
                    * static_cast < PrecisionType > (  vcl_sqrt(  l * ( l + 1 )  )  )
                    * (  pNALPs[ ShIndex( l, -1 ) ]  -  pNALPs[ ShIndex( l, 1 ) ]  )
                    / 2.0;
             

      // the SECOND-order derivative of the OD
      //
      // A1 = sqrt(  l * ( l + 1 ) * ( l - 1 ) * ( l + 2 )  )
      // B1 = Q( l, -2 )
      //
      // A2 = sqrt(  l * ( l + 1 ) * ( l - 1 ) * ( l + 2 )  ) = A1
      // B2 = Q( l,  2 )
      //
      // A3 = l * ( l + 1 ) + l * ( l + 1 ) = 2 * l * ( l + 1 )  
      // B3 = Q( l,  0 )
      //
      // total = ( A1 * B1 + A2 * B2 - A3 * B3 ) / 4
      //       = (  A1 * ( B1 + B2 )  -  A3 * B3  ) / 4
      //
      ( *d2SH_del2 ) += static_cast < PrecisionType > (  shCoefs[ ShIndex( l, 0 ) ]  )
                      * (   static_cast < PrecisionType >
                            (  vcl_sqrt(  l * ( l + 1 ) * ( l - 1 ) * ( l + 2 )  )  )
                          * (   pNALPs[ ShIndex( l, -2 ) ]
                              + pNALPs[ ShIndex( l,  2 ) ]
                            )
                          - 2 * l * ( l + 1 )
                          * pNALPs[ ShIndex( l, 0 ) ]
                        ) / 4.0;
    } // endif: a non-ZERO (even) order l

  } // endfor: each EVEN order l

  
  // Updated below are the three values by accumulating the contributions from the spherical
  // harmonics with non-ZERO phase factors. Similarly updated are the three azimuth-related
  // derivatives: dSH_daz, dSH_deldaz, and d2SH_daz2.
  //
  for ( l =  2; l <= MaxOrder; l += 2 )
  for ( m = -l; m <= l;        m ++   )
  { 
    // skip phase factor 0 since we have addressed it above as a special case
    //
    if ( m == 0 ) continue;

    PrecisionType cosAz = static_cast < PrecisionType > (  vcl_cos ( m * azimuth )  );
    PrecisionType sinAz = static_cast < PrecisionType > (  vcl_sin ( m * azimuth )  );
    
    // the OD value
    tmpVal        = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                  * pNALPs[ ShIndex( l, m ) ];
    tmpVal       *= ( m < 0 ) ? cosAz : sinAz;
    ( *ODValue ) += tmpVal;

    // ----------------------------------------------------------------------
    // the FIRST-order derivative of OD ---- dSH( l, m, theta, phi ) / dTheta
    //
    PrecisionType derive1 = 0.0;
    
    // + ( A1 * B1 )
    if (  ( -l ) <= ( m - 1 )  &&  ( m - 1 ) <= l  )
    {
      derive1 += static_cast < PrecisionType > (  vcl_sqrt(  ( l + m ) * ( l - m + 1 )  )  )
               * pNALPs[ ShIndex( l, m - 1 ) ];
    }
        
    // - ( A2 * B2 )
    if (  ( -l ) <= ( m + 1 )  &&  ( m + 1 ) <= l  )
    {
      derive1 -= static_cast < PrecisionType > (  vcl_sqrt(  ( l - m ) * ( l + m + 1 )  )  )
               * pNALPs[ ShIndex( l, m + 1 ) ];
    }

    derive1      /= -2.0;  // the original sign is NEGATED here xxx -
    tmpVal        = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                  * derive1;
    tmpVal       *= ( m < 0 ) ? cosAz : sinAz;
    ( *dSH_del ) += tmpVal;
    //
    // the FIRST-order derivative of OD ---- dSH( l, m, theta, phi ) / dTheta
    // ----------------------------------------------------------------------


    // -------------------------------------------------------------------------
    // the SECOND-order derivative of OD ---- dSH2( l, m, theta, phi ) / dTheta2
    //
    // - ( A3 * B3 ) always needs to be computed
    PrecisionType derive2 = - (  ( l + m ) * ( l - m + 1 ) + ( l - m ) * ( l + m + 1 )  )
                          * pNALPs[ ShIndex( l, m ) ];
    
    // + ( A1 * B1 )
    if (  ( -l ) <= ( m - 2 )  &&  ( m - 2 ) <= l  )
    {
      derive2 += static_cast < PrecisionType >
                 (  vcl_sqrt(  ( l + m ) * ( l - m + 1 ) * ( l + m - 1 ) * ( l - m + 2 )  )  )
               * pNALPs[ ShIndex( l, m - 2 ) ];
    }
    
    // + ( A2 * B2 )
    if (  ( -l ) <= ( m + 2 )  &&  ( m + 2 ) <= l  )
    {
      derive2 += static_cast < PrecisionType >
                 (  vcl_sqrt(  ( l - m ) * ( l + m + 1 ) * ( l - m - 1 ) * ( l + m + 2 )  )  )
               * pNALPs[ ShIndex( l, m + 2 ) ];
    }

    derive2        /= 4.0;
    tmpVal          = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                    * derive2;
    tmpVal         *= ( m < 0 ) ? cosAz : sinAz;
    ( *d2SH_del2 ) += tmpVal;
    //
    // the SECOND-order derivative of OD ---- dSH2( l, m, theta, phi ) / dTheta2
    // -------------------------------------------------------------------------

      
    // below are the derivatives relative to azimuth
    //
    if ( atPole )
    { 
      tmpVal        = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                    * derive1;
      tmpVal       *= ( m < 0 ) ? ( -sinAz )
                                :    cosAz;  // NO m-scaling here
      ( *dSH_daz ) += tmpVal;
    }
    else
    { 
      // so far, derive1 refers to dSH_del (though a negated version, see the above)
      // and let's add the derivative relative to azimuth
      //
      tmpVal            = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                        * derive1;
      tmpVal           *= ( m < 0 ) ? ( - m * sinAz )      // from cos
                                    : (   m * cosAz );     // from sin
      ( *d2SH_deldaz ) += tmpVal;

      // the value shared by 1st- and 2nd-order pure derivatives relative to azimuth
      //
      PrecisionType comVal   = static_cast < PrecisionType > (  shCoefs[ ShIndex( l, m ) ]  )
                             * pNALPs[ ShIndex( l, m ) ];

      // the FIRST-order derivative of OD relative to azimuth
      //
      tmpVal          = comVal;
      tmpVal         *= ( m < 0 ) ? ( - m * sinAz )        // from cos
                                  : (   m * cosAz );       // from sin
      ( *dSH_daz )   += tmpVal;

      // the SECOND-order derivative of OD relative to azimuth
      //
      tmpVal          = comVal;
      tmpVal         *= ( m < 0 ) ? ( - m * m * cosAz )    // from cos
                                  : ( - m * m * sinAz );   // from sin
      ( *d2SH_daz2 ) += tmpVal;
    }

  } // endfor: double loops --- each ( l, m )


  delete [] pNALPs;
  pNALPs = NULL;

  
  if ( !atPole ) 
  {
    // sinEle = sin( elevation ) = the radius of the xy-circle on which the azimuth angle is based
    // 
    // for any non-pole point, let's make the azimuth-related change independent of the elevation
    // 
    // NOTE: this pre-division by sinEle, sin(elevation), compensates to some degree for the lack of 
    //	     sin(elevation) when determining the dx and dy at the end of the peaking finding function
    //       
    ( *dSH_daz     ) /= sinEle;
    ( *d2SH_deldaz ) /= sinEle;
    ( *d2SH_daz2   ) /= sinEle * sinEle;
  }

} // end function: GetODValueAndDerivatives( ... )


// Given an initial/sample/guess normalized 3D vector (normVec0), this function 
// uses the Newton-Raphson Gradient Descent method to find a peak and returns the
// the normalized peak orientation (normVec1), the associated OD value (ODValue),
// and the status (0: failure; 1: success).
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
int PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::FindPeak( const InputComponentType normVec0[3], 
            OutputComponentType normVec1[3], OutputComponentType * ODValue, const int numIters )
{
  int           beFound = 0;
  PrecisionType tempVal,   tempVec[3];
  PrecisionType azimuth,   elevation;
  PrecisionType dSH_del,   dSH_daz,     derivsMag;
  PrecisionType d2SH_del2, d2SH_deldaz, d2SH_daz2;
  PrecisionType del,       daz;
  PrecisionType dSH_dt,    d2SH_dt2,    dt;

  *ODValue    = 0.0;
  normVec1[0] = static_cast < OutputComponentType > ( normVec0[0] );
  normVec1[1] = static_cast < OutputComponentType > ( normVec0[1] );
  normVec1[2] = static_cast < OutputComponentType > ( normVec0[2] );

  tempVal     = 0.0;
  tempVec[0]  = static_cast < PrecisionType > ( normVec0[0] );
  tempVec[1]  = static_cast < PrecisionType > ( normVec0[1] );
  tempVec[2]  = static_cast < PrecisionType > ( normVec0[2] );
  
  // use numIters iterations in an effort to find a peak
  for ( int i = 0; i < numIters; i ++ )
  {
    // obtain the values of the OD and derivatives based on the current direction
    azimuth   = static_cast < PrecisionType > (  vcl_atan2( tempVec[1], tempVec[0] )  );
    elevation = static_cast < PrecisionType > (  vcl_acos ( tempVec[2] )  );
    GetODValueAndDerivatives( elevation,  azimuth, 
                              &tempVal, 
                              &dSH_del,   &dSH_daz,
                              &d2SH_del2, &d2SH_deldaz, &d2SH_daz2 );

    // normalize the two 1st-order derivatives
    derivsMag = static_cast < PrecisionType >
                (  vcl_sqrt( dSH_del * dSH_del + dSH_daz * dSH_daz )  );
    if ( derivsMag == 0.0 ) break; // epsilon?
    
    del       = dSH_del / derivsMag;
    daz       = dSH_daz / derivsMag;

    // the first-order TOTAL derivative
    dSH_dt    = dSH_del * del + dSH_daz * daz;

    // the second-order TOTAL derivative
    d2SH_dt2  = d2SH_del2 * del * del + d2SH_deldaz * ( 2.0 * del * daz )
              + d2SH_daz2 * daz * daz;
    if ( d2SH_dt2 == 0.0 ) break;  // epsilon ?

    // the combined change in direction over the sphere
    // Newton-Raphson: x[n+1] = x[n] + displacement = x[n] - ( df / d2f )
    dt        = - dSH_dt / d2SH_dt2;
    if ( dt < 0.0 || dt > MAX_DIR_CHANGE ) dt = MAX_DIR_CHANGE;

    // the change in elevation and that in azimuth
    del *= dt;
    daz *= dt;

    // spherical coordinates
    // 
    // x = sin(el) * cos(az)
    // y = sin(el) * sin(az)
    // z = cos(el)
    //
    // dx = [ cos(el) * del ] * cos(az) + sin(el) * [ -sin(az) * daz ] 
    //    = del * cos(az) * cos(el) - daz * sin(az) * sin(el)  ---  inconsistency
    //
    // dy = [ cos(el) * del ] * sin(az) + sin(el) * [  cos(az) * daz ] 
    //    = del * sin(az) * cos(el) + daz * cos(az) * sin(el)  ---  inconsistency
    //
    // dz = -sin(el) * del = -del * sin(el)  ---  consistency
    // 
    // the two inconsistencies are caused by the pre-division, 
    // in GetODValueAndDerivatives( ... ), of sin(elevation) for the 
    // azimuth-related derivatives
    // 
    tempVec[0] += del * static_cast < PrecisionType > (  vcl_cos(  azimuth  )  )
                      * static_cast < PrecisionType > (  vcl_cos( elevation )  )
                - daz * static_cast < PrecisionType > (  vcl_sin(  azimuth  )  );
    tempVec[1] += del * static_cast < PrecisionType > (  vcl_sin(  azimuth  )  )
                      * static_cast < PrecisionType > (  vcl_cos( elevation )  ) 
                + daz * static_cast < PrecisionType > (  vcl_cos(  azimuth  )  );
    tempVec[2] +=-del * static_cast < PrecisionType > (  vcl_sin( elevation )  );

    // normalize the new vector
    PrecisionType  vecMag = static_cast < PrecisionType >
                            (  vcl_sqrt( tempVec[0] * tempVec[0] +
                                         tempVec[1] * tempVec[1] +
                                         tempVec[2] * tempVec[2] )  );
    if ( vecMag == 0.0 ) // epsilon ?
    {
      tempVec[0] = tempVec[1] = tempVec[2] = 0.0;
      break;
    }
    else
    {
      PrecisionType scale = 1.0 / vecMag;
      tempVec[0] *= scale;
      tempVec[1] *= scale;
      tempVec[2] *= scale;
    }

    // convergence detection
    if ( dt < MIN_DIR_CHANGE )  {  beFound = 1;  break;  }
  }

  normVec1[0] = static_cast < OutputComponentType > ( tempVec[0] );
  normVec1[1] = static_cast < OutputComponentType > ( tempVec[1] );
  normVec1[2] = static_cast < OutputComponentType > ( tempVec[2] );
  *ODValue    = static_cast < OutputComponentType > ( tempVal    );
  
  return beFound;
} // end function: FindPeak( ... )


// This function calls FindPeak( ... ) with a set of sample directions to locate 
// all potential/raw peaks, which are sorted (in decreasing order by the OD value)
// and internally stored in an IVAR of type std::vector. Returned is the number
// of the extracted raw peaks.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
int PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::FindRawPeaks()
{
  if ( beReady == 1 ) return sortedPeaks.size() / 4;

  sortedPeaks.clear();

  int                 i, j;
  int                 beFound,    xtracts;
  double              theta, phi;
  PrecisionType       sinTheta,   odValue;
  PrecisionType       vector0[3], vector1[3];
  OutputComponentType tempVec[3];

  
  for ( j = 0; j < NUM_SAMPLEDIRS; j ++)
  { 
    // get the sample direction as a 3D vector
    int  index = ( j << 1 );
    phi        = double(  SAMPLE_DIRS[ index     ]  );
    theta      = double(  SAMPLE_DIRS[ index + 1 ]  );
    
    sinTheta   = static_cast < PrecisionType > (  vcl_sin( theta )  );
    vector0[0] = static_cast < PrecisionType >
                 (    sinTheta * static_cast < PrecisionType > 
                                 (  vcl_cos( phi )  )    );
    vector0[1] = static_cast < PrecisionType >
                 (    sinTheta * static_cast < PrecisionType > 
                                 (  vcl_sin( phi )  )    );
    vector0[2] = static_cast < PrecisionType > (  vcl_cos( theta )  );

    /*/ normalization (for more accuracy?)
    PrecisionType scale = 1.0 / static_cast < PrecisionType >
                                (  vcl_sqrt( vector0[0] * vector0[0] +
                                             vector0[1] * vector0[1] +
                                             vector0[2] * vector0[2] )  );
    vector0[0] *= scale;
    vector0[1] *= scale;
    vector0[2] *= scale; //*/

    // find a peak, if any, based on (closest to) the sample direction
    if (  FindPeak( vector0, vector1, &odValue, 50 )  ==  0  )  continue;

    // check if this peak is highly close to an already extracted one
    beFound = 0;
    xtracts = sortedPeaks.size() / 4; // 4 components per peak: vx, vy, vz, mag
    for (  i = 0;  ( i < xtracts ) && ( beFound == 0 );  i ++  )
    {
      // get the peak as a 3D vector
      index      = ( i << 2 );
      tempVec[0] = sortedPeaks[ index     ];
      tempVec[1] = sortedPeaks[ index + 1 ];
      tempVec[2] = sortedPeaks[ index + 2 ];

      // check how close this existing peak is to the new
      // one so as to discard duplicate and paired peaks  
      // 
      PrecisionType vecsDot = tempVec[0] * vector1[0] + 
                              tempVec[1] * vector1[1] + 
                              tempVec[2] * vector1[2];
      if ( vecsDot >=  MAX_DOTPRODUCT || 
           vecsDot <= -MAX_DOTPRODUCT ) beFound = 1;
    }

    // register the new peak if it differs sufficiently from the existing ones
    if ( beFound == 0 )
    {
      sortedPeaks.push_back(  static_cast < OutputComponentType > ( vector1[0] )  );
      sortedPeaks.push_back(  static_cast < OutputComponentType > ( vector1[1] )  );
      sortedPeaks.push_back(  static_cast < OutputComponentType > ( vector1[2] )  );
      sortedPeaks.push_back(  static_cast < OutputComponentType > ( odValue    )  );
    }
      
  } // endfor: each sample direction


  // sort the extracted peaks by the OD value in decreasing order
  xtracts = sortedPeaks.size() / 4;
  for ( j = 0; j < xtracts - 1; j ++ )
  {
    int   index0 = ( j << 2 );
    for ( i = j + 1; i < xtracts; i ++ )
    {
      int index1 = ( i << 2 );
      if (  sortedPeaks[ index0 + 3 ]  <  sortedPeaks[ index1 + 3 ]  )
      {
        for ( int k = 0; k < 4; k ++ )
        {
          OutputComponentType tempVal = sortedPeaks[ index0 + k ];
          sortedPeaks[ index0 + k ]   = sortedPeaks[ index1 + k ];
          sortedPeaks[ index1 + k ]   = tempVal;
        }
      }
    }
  }


  // remove fake peaks --- those with non-positive OD values
  if ( toCheck == 1 && xtracts > 0 ) // ON except for debugging
  {     
    // locate the first nonpositive-OD peak
    for ( j = 0; j < xtracts; j ++ )
    {
      if (    sortedPeaks[  ( j << 2 ) + 3  ]    <=    0.0    )  break;
    }
    
    // now j indicates the number of POSITIVE UNI-DIRECTIONAL peaks
    // stored in the contiguous former part
    sortedPeaks.resize( j << 2 );       // four components per peak
  }
  

  beReady = 1;
  return sortedPeaks.size() / 4;
} // end function: FindRawPeaks()


// Get the raw peaks, each with a normalized 3D direction (vx, vy, vz) and the
// associated OD value --- a 4-component tuple.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
int PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::GetRawPeaks( OutputType & peaks )
{
  peaks.clear();
  int numRaws = FindRawPeaks();
  int numVals = ( numRaws << 2 );

  for ( int i = 0; i < numVals; i ++ )
  {
    peaks.push_back( sortedPeaks[i] );
  }

  return numRaws;
}


// Get the range of the OD values at raw peaks (max followed by min). Returned is
// a flag: 0 for failure (no any peak at all) and non-zero (number of raw peaks)
// for success.
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
int PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::GetRangeOfRawPeaks( OutputComponentType maxMin[2] )
{
  maxMin[0]     = 
  maxMin[1]     =-1000.0;

  int  numPeaks = FindRawPeaks();
  if ( numPeaks > 0 )
  {
    maxMin[0] = sortedPeaks[3];
    maxMin[1] = sortedPeaks[  (  ( numPeaks - 1 ) << 2  )  +  3  ];
  }

  return numPeaks;
}


// Get the normalized 3D direction (normVec) and the associated OD value (odValue) 
// of the peakIdx-th (0-based) largest peak, either using the range limit (1) or 
// not (0). Returned is a flag indicating whether the peak is available (peakIdx 
// is valid: 1) or not (0).
template < typename TPixelType, typename TOutputComponentType, typename TPrecisionType >
int PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
::GetPeak( int peakIdx, OutputComponentType normVec[3], 
           OutputComponentType & odValue, int useRange )
{ 
  // init the result --- necessary for the zero-padding case 
  normVec[0] = 0.0;
  normVec[1] = 0.0;
  normVec[2] = 0.0;
  odValue    = 0.0;

  int  numPeaks = FindRawPeaks();
  if ( numPeaks < 1 ) return 0;

  int  peakIdx0 = 0;            // the most left  index
  int  peakIdx1 = numPeaks - 1; // the most right index

  if ( useRange )
  {
    int   i;
    for ( i = 0; i < numPeaks; i ++ )      // max: the most left
    { 
      if (  rangeOfInterest[0]  >=  sortedPeaks[ ( i << 2 ) + 3 ]  )  break;
    }
    peakIdx0 = i;

    for ( i = numPeaks - 1; i >= 0; i -- ) // min: the most right
    {
      if (  rangeOfInterest[1]  <=  sortedPeaks[ ( i << 2 ) + 3 ]  )  break;
    }
    peakIdx1 = i;
  }

  if ( peakIdx0 >= numPeaks || peakIdx1 < 0 || 
       peakIdx0 >  peakIdx1 || peakIdx0 + peakIdx > peakIdx1 ) return 0;

  int  theIndex = ( peakIdx0 + peakIdx ) << 2;
  normVec[0] = sortedPeaks[ theIndex     ];
  normVec[1] = sortedPeaks[ theIndex + 1 ];
  normVec[2] = sortedPeaks[ theIndex + 2 ];
  odValue    = sortedPeaks[ theIndex + 3 ];

  return 1;
} // end function: GetPeak( ... )


} // end namespace


#endif //__itkPeakFindingCalculator_hxx
