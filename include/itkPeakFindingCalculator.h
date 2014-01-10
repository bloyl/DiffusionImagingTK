/**
 * @file  itkPeakFindingCalculator.h
 * @brief -
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkPeakFindingCalculator_h
#define __itkPeakFindingCalculator_h


// require access to itkSymRealSphericalHarmonicRep::K( int l, int m )
#include <itkSymRealSphericalHarmonicRep.h>


namespace itk
{

template < typename TPixelType, 
           typename TOutputComponentType = double, 
           typename TPrecisionType       = double >
class ITK_EXPORT PeakFindingCalculator: public Object
{
public:

  // Standard class typedefs.
  typedef PeakFindingCalculator< TPixelType, TOutputComponentType, TPrecisionType >
                                                        Self;
  typedef SmartPointer < Self >                         Pointer;
  typedef SmartPointer < const Self >                   ConstPointer;

  // Input and output types
  // 
  // Input:  a total of ( MaxOrder + 1 ) * ( MaxOrder + 2 ) / 2 spherical harmonic
  //         coefficients as a representation of an OD function (either FOD or ODF).
  //
  // Output: the peaks (local extrema) of the OD function in the form of normalized
  //         3D vectors/directions and the associated magnitudes (OD values).
  //
  // NOTE:   Each fiber results in a pair of mutually-opposite directions, both with
  //         positive OD values. This rule is used to remove fake/artifact peaks.
  //
  typedef TPixelType                                    PixelType;
  typedef TPrecisionType                                PrecisionType;
  typedef TOutputComponentType                          OutputComponentType;
  typedef typename PixelType::ComponentType             InputComponentType;

  itkStaticConstMacro( MaxOrder,  unsigned int, PixelType::MaxOrder  );
  itkStaticConstMacro( Dimension, unsigned int, PixelType::Dimension );

  typedef itk::Vector < InputComponentType, Dimension > InputType;
  typedef std::vector < OutputComponentType >           OutputType;

  typedef  struct
  { 
    TPrecisionType   weight0;
    TPrecisionType * lutPtr0;
    TPrecisionType   weight1;
    TPrecisionType * lutPtr1;
  } Interp;
  
  // Method for creation through the object factory.
  itkNewMacro(Self);  

  // Run-time type information (and related methods).
  itkTypeMacro( PeakFindingCalculator, Object );


  // ===============================================================================


  // Turn on /off the use of an LUT for accelrated NALP computation
  void SetAcceleration ( int lutUsed = 1 )  {  LUTused = lutUsed;  }

  // Turn on/off the peak filter (on by default), which is used to remove fake peaks
  // -- those with non-positive OD values
  void SetPeakFiltering( int tocheck = 1 )  {  toCheck = tocheck;  }

  // Set a spherical harmonic coefficient specified by an index [ 0, Dimension - 1 ].
  // No range check is performed here and the caller needs to guarantee a valid index.
  inline
  void SetCoefficient( int index, InputComponentType value )
       {  this->shCoefs[ index ] = value;  beReady = 0;  }


  // Set a specified number (Dimension) of spherical harmonic coefficients.
  void SetCoefficients( InputType & coefs )  
       {  this->shCoefs = coefs;  beReady = 0;  }

  // Set a range (max followed by min) of interest (in terms of the OD value) 
  // within which to access the peaks.
  void SetRangeOfInterest( OutputComponentType maxMin[2] );

  // Get the number of raw peaks.
  int  GetNumberOfRawPeaks()  {  return FindRawPeaks();  }

  // Get the range of the OD values at raw peaks (max followed by min). Returned is
  // a flag: 0 for failure (no any peak at all) and non-zero (number of raw peaks)
  // for success.
  int  GetRangeOfRawPeaks( OutputComponentType maxMin[2] );

  // Given an initial/sample/guess normalized 3D vector (normVec0), this function 
  // uses the Newton-Raphson Gradient Descent method to find a peak and returns the
  // the normalized peak orientation (normVec1), the associated OD value (ODValue),
  // and the status (0: failure; 1: success).
  int  FindPeak( const InputComponentType normVec0[3], 
                 OutputComponentType normVec1[3], OutputComponentType * ODValue, const int numIter );

  // Get the raw peaks, each with a normalized 3D direction (vx, vy, vz) and the
  // associated OD value --- a 4-component tuple.
  int  GetRawPeaks( OutputType & peaks );

  // Get the normalized 3D direction (normVec) and the associated OD value (odValue) 
  // of the peakIdx-th (0-based) largest peak, either using the range limit (1) or 
  // not (0). Returned is a flag indicating whether the peak is available (peakIdx 
  // is valid: 1) or not (0). Upon unavailability, normVec = (0.0, 0.0, 0.0) and
  // odValue = 0.0 are returned for zero-padding purposes.
  int  GetPeak( int peakIdx, OutputComponentType normVec[3], 
                OutputComponentType & odValue, int useRange = 0 );


protected: 

  PeakFindingCalculator();
  ~PeakFindingCalculator();


private:

  PeakFindingCalculator( const Self & ); //purposely not implemented
  void operator = ( const Self & );      //purposely not implemented

  // generate an LUT for accelerated computation of NALP values
  void GenerateNalpLUT();

  // set up the context for interpolating NALP values across elevations
  void SetInterpolator( TPrecisionType eleVal, Interp * interp );

  // get an NALP value through elevation-based linear interpolation
  inline TPrecisionType NALP( const Interp * interp, int l, int m ) 
  { 
    int            shIndex = ShIndex( l, m );
    TPrecisionType nalpVal = interp->lutPtr0[ shIndex ] * interp->weight0
                           + interp->lutPtr1[ shIndex ] * interp->weight1;
    return nalpVal;
  }

  // Given EVEN order l ( l >= 0 ) and phase factor m ( -l <= m <= l ), return
  // the ZERO-based index of a spherical harmonic function / coefficient in the
  // representation scheme.
  inline unsigned int ShIndex( int l, int m )
  {  return  (    (  l * ( l + 1 )  )  >>  1    )    +    m;  }

  // The Normalization factor for spherical harmonic (l, m). xxx
  PrecisionType K( int l, int m );

  // The Normalized Associated Legender Polynomial (NALP) value for spherical 
  // harmonic (l, m, theta).
  PrecisionType NALP( int l, int m, PrecisionType elevation );

  // This function computes the OD as well as 1st-order and 2nd-order derivatives
  // along a specific direction (elevation, azimuth).
  // 
  // As opposed to Tournier's representation, cos is assigned to negative phases
  // and sin to positive phases in class itkSymRealSphericalHarmonicRep, which
  // results in a different set of equations for calculating derivatives and OD
  // values.
  void GetODValueAndDerivatives
       ( PrecisionType   elevation,   PrecisionType   azimuth,
         PrecisionType * ODValue,     PrecisionType * dSH_del,
         PrecisionType * dSH_daz,     PrecisionType * d2SH_del2,
         PrecisionType * d2SH_deldaz, PrecisionType * d2SH_daz2 );

  // This function calls FindPeak( ... ) with a set of sample directions to locate 
  // all potential/raw peaks, which are sorted (in decreasing order by the OD value)
  // and internally stored in an IVAR of type std::vector. Returned is the number
  // of the extracted raw peaks.
  int  FindRawPeaks();

  int                 beReady;            // are raw peaks available?
  int                 toCheck;            // to remove fake peaks?
  int                 LUTused;            // is an LUT used for accelerated NALP computation?
  InputType           shCoefs;            // spherical harmonic coefficients
  OutputType          sortedPeaks;        // peaks sorted by the OD amplitude
  TPrecisionType   *  nalpLUT;            // an LUT for accelerated NALPs computation
  OutputComponentType rangeOfInterest[2]; // NOTE: max followed by min

};
  

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPeakFindingCalculator.hxx"
#endif

  
#endif
