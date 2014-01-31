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
