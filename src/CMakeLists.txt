

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(VNL_MATRIX_TEMPLATES)
foreach(n IN ITEMS 15 28 45 66 91 120 153 190 231)
  list(APPEND VNL_MATRIX_TEMPLATES "vtkMatrixTemplates/vnl_matrix_fixed.${n}.${n}.cxx")
endforeach()

# add_library(DWICore
#   ${VNL_MATRIX_TEMPLATES}
#   itkReplaceSpecialFunctions.cxx
#   # sbiaBasicUtilities.cxx
#   # itkNumericTraitsVectorPixel2.cxx
#   )
# target_link_libraries(DWICore ${ITK_LIBRARIES})

set(${itk-module}_SRC
  ${VNL_MATRIX_TEMPLATES}
  itkReplaceSpecialFunctions.cxx
)

add_library(${itk-module} ${${itk-module}_SRC})
target_link_libraries(${itk-module}  ${${itk-module}_LIBRARIES})
itk_module_target(${itk-module})
