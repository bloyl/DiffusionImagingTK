itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}UnitTests
  itkSymRealSphericalHarmonicRepTest.cxx
  itkPixelReorientationOperatorTest.cxx
  itkDiffusionModelCalculatorTest.cxx
  itkPeakFindingCalculatorGridTest.cxx
  itkPeakFindingCalculatorTest.cxx
  itkTrackerDirectionPickerDWITest.cxx
  itkTrackerDirectionPickerRSHTest.cxx
  itkTrackerDirectionPickerDTITest.cxx
  # itkFiberGeneratorTest.cxx
  # itkFiberTrackingManagerTest.cxx
#  itkFiberTrackingManagerTest2.cxx
#  itkTestBoundingBox.cxx 
)

CreateTestDriver(${itk-module}Unit  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}UnitTests}")

set(TestNamePrefix "${itk-module}_")

itk_add_test(NAME "${TestNamePrefix}itkSymRealSphericalHarmonicRepTest"
  COMMAND ${itk-module}UnitTestDriver itkSymRealSphericalHarmonicRepTest )

itk_add_test(NAME "${TestNamePrefix}itkPixelReorientationOperatorTest"
  COMMAND ${itk-module}UnitTestDriver itkPixelReorientationOperatorTest )

itk_add_test(NAME "${TestNamePrefix}tkDiffusionModelCalculatorTest"
  COMMAND ${itk-module}UnitTestDriver itkDiffusionModelCalculatorTest )

itk_add_test(NAME "${TestNamePrefix}itkPeakFindingCalculatorGridTest"
  COMMAND ${itk-module}UnitTestDriver itkPeakFindingCalculatorGridTest )

itk_add_test(NAME "${TestNamePrefix}itkPeakFindingCalculatorTest"
  COMMAND ${itk-module}UnitTestDriver itkPeakFindingCalculatorTest )

itk_add_test(NAME "${TestNamePrefix}itkTrackerDirectionPickerDWITest"
  COMMAND ${itk-module}UnitTestDriver itkTrackerDirectionPickerDWITest )

itk_add_test(NAME "${TestNamePrefix}itkTrackerDirectionPickerDTITest"
  COMMAND ${itk-module}UnitTestDriver itkTrackerDirectionPickerDTITest )

itk_add_test(NAME "${TestNamePrefix}itkTrackerDirectionPickerRSHTest"
  COMMAND ${itk-module}UnitTestDriver itkTrackerDirectionPickerRSHTest )

# ##########
# # list of test implementations
# set (
#   DWICoreSOURCES
#     itkPixelReorientationOperatorTest.cxx
#     itkDiffusionModelCalculatorTest.cxx
#     itkPeakFindingCalculatorGridTest.cxx
#     itkSymRealSphericalHarmonicRepTest.cxx
#     itkPeakFindingCalculatorTest.cxx
# ) 

# #     itkVariableLengthVectorCastImageFilterTest.cxx
# #     itkFiberImpulseResponseImageCalculatorTest.cxx
# #     itkFodCsdReconImageFilterTest.cxx
# #     itkOdfReconImageFilterTest.cxx
# #     itkSymRSHPowerImageFilterTest.cxx
# #     itkSymRshReconImageFilterTest.cxx
# # )

# set (LIBS DWICore ${ITK_LIBRARIES})
# if ("${ITK_VERSION_MAJOR}" LESS 4)
#   list (APPEND LIBS ITKNumerics ITKCommon ITKIO)
# endif ()

# # build test driver
# create_test_sourcelist (DWICoreTest_SOURCES DWICoreSourceTestsDriver.cxx ${DWICoreSOURCES})

# # build test driver
# add_executable (DWICoreTest ${DWICoreTest_SOURCES})
# target_link_libraries (DWICoreTest ${LIBS})

# # add tests
# foreach(source IN LISTS DWICoreSOURCES)
#   get_filename_component(test_name ${source} NAME_WE)
#   add_test( ${test_name} "${EXECUTABLE_OUTPUT_PATH}/DWICoreTest" "${test_name}" )
# endforeach()

# ##############################################################################################################
# set (FiberTrackingSOURCES
#   itkTrackerDirectionPickerDWITest.cxx
#   itkTrackerDirectionPickerRSHTest.cxx
#   itkTrackerDirectionPickerDTITest.cxx
#   itkFiberGeneratorTest.cxx
#   itkFiberTrackingManagerTest.cxx
# #  itkFiberTrackingManagerTest2.cxx
#   itkTestBoundingBox.cxx
#   )

# set(LIBS DWICore; ${ITK_LIBRARIES}; ${VTK_LIBRARIES}; )

# if( "${ITK_VERSION_MAJOR}" LESS 4 )
#   list(APPEND LIBS ITKNumerics ITKCommon ITKIO)
# endif( "${ITK_VERSION_MAJOR}" LESS 4 )

# # build test driver
# create_test_sourcelist (FiberTrackingTest_SOURCES FiberTrackingTestsDriver.cxx ${FiberTrackingSOURCES})

# # build test driver
# add_executable (FiberTrackingTests ${FiberTrackingTest_SOURCES})
# target_link_libraries (FiberTrackingTests ${LIBS})

# # add tests
# foreach(source IN LISTS FiberTrackingSOURCES)
#   get_filename_component(test_name ${source} NAME_WE)
#   add_test( ${test_name} "${EXECUTABLE_OUTPUT_PATH}/FiberTrackingTests" "${test_name}" "${TESTING_INPUT_DATA_PATH}" "${TESTING_OUTPUT_DIR}")
# endforeach()
