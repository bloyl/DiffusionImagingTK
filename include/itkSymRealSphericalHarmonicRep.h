/**
 * @file  itkSymRealSphericalHarmonicRep.h
 * @brief Supplies the basic representation for real valued antipodally symmetric spherical functions.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#ifndef __itkSymRealSphericalHarmonicRep_h
#define __itkSymRealSphericalHarmonicRep_h

// Undefine an eventual SymRealSphericalHarmonicRep macro
#ifdef SymRealSphericalHarmonicRep
#  undef SymRealSphericalHarmonicRep
#endif

#include <itkIndent.h>
#include <itkFixedArray.h>
#include <itkMatrix.h>
#include <itkNumericTraits.h>
#include <itkArray2D.h>
#include <itkSymmetricEigenAnalysis.h>
#include <itkVectorContainer.h>
#include <vnl/vnl_vector.h>
#include <vnl/vnl_vector_fixed.h>

#include <itkReplaceSpecialFunctions.h>

namespace itk
{

template
< typename TComponent,
  unsigned int TMaxOrder=4
//  ,typename TBasisType= itk::RealSymSphericalHarmonicBasis< TMaxOrder >
>
class SymRealSphericalHarmonicRep: public
      FixedArray<TComponent,(TMaxOrder+1)*(TMaxOrder+2)/2>
{
public:

  /** Standard class typedefs. */
  typedef SymRealSphericalHarmonicRep  Self;
  typedef FixedArray<TComponent,(TMaxOrder+1)*(TMaxOrder+2)/2> Superclass;

  /** Runtime information support. */
  itkTypeMacro(SymRealSphericalHarmonicRep, FixedArray);

  /** Dimension of the vector space. */
  itkStaticConstMacro(MaxOrder, unsigned int, TMaxOrder);
  itkStaticConstMacro(Dimension, unsigned int, (TMaxOrder+1)*(TMaxOrder+2)/2);

  /**Dimension of unique Orders since Max Order must be even. */
  itkStaticConstMacro(NumberOfOrders, unsigned int, TMaxOrder / 2 + 1);

  typedef vnl_vector_fixed<int,2>                     LmVector;

  /** Convenience typedefs. */
  typedef FixedArray<TComponent, itkGetStaticConstMacro(Dimension)> BaseArray;

  /**  Define the component type. */
  typedef TComponent ComponentType;
  typedef typename Superclass::ValueType ValueType;
  typedef typename NumericTraits<ValueType>::RealType AccumulateValueType;
  typedef typename NumericTraits<ValueType>::RealType RealValueType;

  /**  Define the Gradient Direction Type. */
  typedef vnl_vector_fixed< double, 3 >               GradientDirectionType;

  /** Container to hold gradient directions of the 'n' DW measurements */
  typedef VectorContainer< unsigned int,
          GradientDirectionType >                     GradientDirectionContainerType;

  typedef vnl_matrix<double>                          RshBasisMatrixType;

  /** Default constructor. */
  SymRealSphericalHarmonicRep()
    {
    if ( (MaxOrder % 2) != 0 )
      {
      itkGenericExceptionMacro( << "Symetric real spherical harmonic representations are only of even order");
      }

    this->Fill(0);
    }

  SymRealSphericalHarmonicRep (const ComponentType& r)
    {
    if ( (MaxOrder % 2) != 0 )
      {
      itkGenericExceptionMacro( << "Symetric real spherical harmonic representations are only of even order!");
      }
    this->Fill(r);
    }

  typedef ComponentType ComponentArrayType[ itkGetStaticConstMacro(Dimension) ];

  /** Pass-through constructor for the Array base class. */
  SymRealSphericalHarmonicRep(const Self& r): BaseArray(r)
    {
    if ( (MaxOrder % 2) != 0 )
      {
      itkGenericExceptionMacro( << "Symetric real spherical harmonic representations are only of even order!");
      }
    }

  SymRealSphericalHarmonicRep(const ComponentArrayType r): BaseArray(r)
    {
    if ( (MaxOrder % 2) != 0 )
      {
      itkGenericExceptionMacro( << "Symetric real spherical harmonic representations are only of even order!");
      }
    }

  /** Templated constructor */
  template < typename TCoordRepB >
  SymRealSphericalHarmonicRep( const SymRealSphericalHarmonicRep<TCoordRepB,TMaxOrder> & pa ):
    BaseArray(pa) { };


  /** Pass-through assignment operator for the Array base class. */
  Self& operator= (const Self& r);
  Self& operator= (const ComponentType& r);
  Self& operator= (const ComponentArrayType r);

  /**
   * Assigment from a vnl_vector
   */
  template <typename VectorComponentType>
  Self& operator= (const vnl_vector<VectorComponentType> r )
  {
    ///check the size...
    if (r.size() != Dimension)
    {
      itkGenericExceptionMacro( << "Assignment must assign every element of array:"
          << "Expecting " << Dimension << "elements. But only received " << r.size()
                     );
    }
    for (unsigned int i=0;i<r.size();i++){
      (*this)[i] = static_cast<ComponentType>(r.get(i));
    }
    return *this;
  }


  /** Aritmetic operations between pixels. Return a new SymRealSphericalHarmonicRep. */
  Self operator+(const Self &vec) const;
  Self operator-(const Self &vec) const;
  const Self & operator+=(const Self &vec);
  const Self & operator-=(const Self &vec);

  /** Arithmetic operations between tensors and scalars */
  Self operator*(const RealValueType & scalar ) const;
  Self operator/(const RealValueType & scalar ) const;
  const Self & operator*=(const RealValueType & scalar );
  const Self & operator/=(const RealValueType & scalar );

  /** Return the number of components. */
  static unsigned int GetNumberOfComponents()
    {
    return itkGetStaticConstMacro(Dimension);
    }

  /** Return the number of Orders. */
  static unsigned int GetNumberOfOrder()
    {
    return itkGetStaticConstMacro(Order);
    }

  /** Return the value for the Nth component. */
  ComponentType GetNthComponent(int c) const
  {
    return this->operator[](c);
  }

  /** Set the Nth component to v. */
  void SetNthComponent(int c, const ComponentType& v)
  {
    this->operator[](c) = v;
  }

  /** Return the value for the l and mth component. */
  ComponentType GetLthMthComponent(int l,int m) const
  {
    if ( (l % 2) != 0 || l > 20 )
    {
      itkGenericExceptionMacro( << "Attempting to extract a componant with an illegal order (l)");
    }
    if ( vcl_abs(m) > l)
    {
      itkGenericExceptionMacro( << "Attempting to Set a componant with an |m| > l");      
    }
    int c = GetJ(l,m);
    return this->operator[](c-1);
  }

  /** Set the Nth component to v. */
  void SetLthMthComponent(int l,int m, const ComponentType& v)
  {
    if ( (l % 2) != 0 || l > 20 )
    {
      itkGenericExceptionMacro( << "Attempting to Set a componant with an illegal order (l)");
    }
    if ( vcl_abs(m) > l)
    {
      itkGenericExceptionMacro( << "Attempting to Set a componant with an |m| > l");      
    }
    int c = GetJ(l,m);
    this->operator[](c-1) = v;
  }
 
  static const RshBasisMatrixType ComputeRshBasis( const GradientDirectionContainerType* );

  BaseArray GetFixedArray() { return (*this); }

  static const LmVector GetLM(unsigned int);

  const RealValueType Evaluate(RealValueType theta, RealValueType phi) const;
  const RealValueType Evaluate(GradientDirectionType Gradient) const;

  const RealValueType ComputeCurvatures(RealValueType theta, RealValueType phi,
            RealValueType& gaussCurvature, RealValueType& meanCurvature, RealValueType& k1, RealValueType& k2) const;
  const RealValueType ComputeCurvatures(GradientDirectionType Gradient,
            RealValueType& gaussCurvature, RealValueType& meanCurvature, RealValueType& k1, RealValueType& k2) const;

  /// Evaluate the jth Basis Function
  static const double Y( int j, double theta, double phi );
  static const double Y( int j, GradientDirectionType Gradient)
  {
    double theta = acos(Gradient[2]);
    double phi   = atan2(Gradient[1],Gradient[0]); // atan2(y,x) = atan(y/x);
    return Y(j,theta,phi);
  }
  
  static double Y( int l, int m, double theta, double phi )
  {
    return Y(GetJ(l,m),theta,phi);
  }

  // Returns the normalization constant for the SH basis function with parameters (l,m).
  // called by itkPeakFindingCalculator
  static const double K( int l, int m );

  static const unsigned int GetJ(int,int);

  void Normalize();
protected:

private:
  /// Returns the normalization constant for the SH basis function with parameters (l,m).
  //static const double K( int l, int m );

};

} // end namespace itk

#include "itkNumericTraitsSymRshPixel.h"
#include "itkSymRealSphericalHarmonicRep.hxx"

#endif
