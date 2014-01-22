
#ifndef __testingUtils_h
#define __testingUtils_h

namespace DiffusionImagingTK_testing
{

template <typename GradientDirectionContainerType>
void loadGradDirs( typename GradientDirectionContainerType::Pointer gradDirs, std::string bvecFile, std::string bvalFile );

template < class DwiImageType4D, class DwiImageType >
void convertDWIimage(typename DwiImageType4D::Pointer inp4d, typename DwiImageType::Pointer outputIm );

}//End Namespace


#endif

#include "testingUtils.hxx"
