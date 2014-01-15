# the top-level README is used for describing this module, just
# re-used it for documentation here
get_filename_component( MY_CURENT_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file( READ "${MY_CURENT_DIR}/README.md" DOCUMENTATION )

# ITK version 4.5 changed it from EXCLUDE_FROM_ALL to EXCLUDE_FROM_DEFAULT
set( _EXCLUDE "EXCLUDE_FROM_ALL" )
if (NOT "${ITK_VERSION_MAJOR}.${ITK_VERSION_MINOR}.${ITK_VERSION_MINOR_PATCH}" VERSION_LESS "4.5")
  set( _EXCLUDE "EXCLUDE_FROM_DEFAULT" )
endif()

list(APPEND ExternalData_URL_TEMPLATES "http://dl.dropboxusercontent.com/u/251948538/%(algo)/%(hash)" )

# itk_module() defines the module dependencies in DiffusionImagingTK
# DiffusionImagingTK depends on ITKCommon
# The testing module in DiffusionImagingTK depends on ITKTestKernel
# and ITKMetaIO(besides DiffusionImagingTK and ITKCore)
# By convention those modules outside of ITK are not prefixed with
# ITK.

# define the dependencies of the include module and the tests
itk_module(DiffusionImagingTK
  DEPENDS
    ITKCommon
    ITKMesh
    ITKImageFunction
    ITKSpatialObjects
    ITKThresholding
  TEST_DEPENDS
    ITKTestKernel
    ITKMetaIO
    ITKIOSpatialObjects
    ITKMesh ## this is needed so tests are correctly linked. Not sure if it should need to be in both places
    ITKSpatialObjects
  DESCRIPTION
    "${DOCUMENTATION}"
  ${_EXCLUDE}
)
