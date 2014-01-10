/**
* @file  itkTrackerDirectionPickerImageBase.hxx
* @brief This class is an abstract class that serves as the base for classes used
* to choose 3 dimensional directions from diffusion model images (ie DTI or FOD images)
*
* Copyright (c) <year> University of Pennsylvania. All rights reserved.
* See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
*
* Contact: SBIA Group <sbia-software at uphs.upenn.edu>
*/

#ifndef __itkTrackerDirectionPickerImageBase_hxx
#define __itkTrackerDirectionPickerImageBase_hxx

#include <itkTrackerDirectionPickerImageBase.h>


namespace itk
{

template <class TInterpType >
TrackerDirectionPickerImageBase<TInterpType>
::TrackerDirectionPickerImageBase()
{
  this->m_Interpolator = NULL;
}

template <class TInterpType >
void
TrackerDirectionPickerImageBase<TInterpType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
  if (this->IsInitialized())
    {
    os << indent << "Is Initialized: " << std::endl;
    os << indent << "Interpolator: " << this->m_Interpolator << std::endl;
    }
  else
    {
    os << indent << "Is *NOT* Initialized: " << std::endl;
    }

  os << std::endl;
}


template <class TInterpType >
bool
TrackerDirectionPickerImageBase<TInterpType>
::IsInitialized() const
{
  if (this->m_Interpolator && this->m_Interpolator->GetInputImage())
    return true;
  else
    return false;
}

template <class TInterpType >
bool
TrackerDirectionPickerImageBase<TInterpType>
::CanEvaluate(PointType pt) const
{
  return (this->m_Interpolator->IsInsideBuffer(pt));
}

//These shouldn't be needed maybe?
template <class TInterpType >
typename TrackerDirectionPickerImageBase<TInterpType>::DirectionType
TrackerDirectionPickerImageBase<TInterpType>
::PickStartingDirection(PointType ) const
{
  DirectionType dir;
  dir[0] = 0;
  dir[1] = 0;
  dir[2] = 0;
  return dir;
}

template <class TInterpType >
typename TrackerDirectionPickerImageBase<TInterpType>::DirectionType
TrackerDirectionPickerImageBase<TInterpType>
::PickNextDirection(PointType , DirectionType ) const
{
  DirectionType dir;
  dir[0] = 0;
  dir[1] = 0;
  dir[2] = 0;
  return dir;
}


}// end namespace itk

#endif
