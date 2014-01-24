itk_module_test()

#${itk-module} will be the name of this module and will not need to be
#changed when this module is renamed.

set(${itk-module}UnitTests
  testingUtils.cxx
  itkTrackerDirectionPickerImageBaseTest.cxx
  itkReplaceSpecialFunctionsTest.cxx
  itkSymRealSphericalHarmonicRepTest.cxx
  itkPixelReorientationOperatorTest.cxx
  itkDiffusionModelCalculatorTest.cxx
  itkPeakFindingCalculatorGridTest.cxx
  itkPeakFindingCalculatorTest.cxx
  itkTrackerDirectionPickerDWITest.cxx
  itkTrackerDirectionPickerRSHTest.cxx
  itkTrackerDirectionPickerDTITest.cxx
  itkFiberGeneratorTest.cxx
  itkFiberTrackingManagerTest.cxx
  # itkFiberTrackingManagerTest2.cxx
  itkDTIReconImageFilterTest.cxx
  itkRSHReconImageFilterTest.cxx
)

CreateTestDriver(${itk-module}Unit  "${${itk-module}-Test_LIBRARIES}" "${${itk-module}UnitTests}")

set(TestNamePrefix "${itk-module}_")

itk_add_test(NAME "${TestNamePrefix}itkTrackerDirectionPickerImageBaseTest"
  COMMAND ${itk-module}UnitTestDriver itkTrackerDirectionPickerImageBaseTest )

itk_add_test(NAME "${TestNamePrefix}itkReplaceSpecialFunctionsTest"
  COMMAND ${itk-module}UnitTestDriver itkReplaceSpecialFunctionsTest )

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

itk_add_test(NAME "${TestNamePrefix}itkFiberGeneratorTest"
  COMMAND ${itk-module}UnitTestDriver itkFiberGeneratorTest DATA{Baseline/sim_dti_orig.nii.gz} )

itk_add_test(NAME "${TestNamePrefix}itkFiberTrackingManagerTest"
  COMMAND ${itk-module}UnitTestDriver itkFiberTrackingManagerTest 
    DATA{Baseline/DWIS_dti-scheme_SNR-30_DTI.nii.gz}
    DATA{Baseline/DWIS_dti-scheme_SNR-30_DTI_FA_mask.nii.gz}
    DATA{Baseline/seeding_regions.nii.gz}
    DATA{Baseline/DWIS_dti-scheme_seedTest_mask.nii.gz}
    DATA{Baseline/DWIS_dti-scheme_singleSeedTest_mask.nii.gz}
    "sim_dti_orig-fibers.meta"
  )

itk_add_test(NAME "${TestNamePrefix}itkDTIReconImageFilterTest"
  COMMAND ${itk-module}UnitTestDriver itkDTIReconImageFilterTest
    DATA{Baseline/template_dwi.nii.gz}
    DATA{Baseline/template_dwi.bvec}
    DATA{Baseline/template_dwi.bval}
    DATA{Baseline/template_dwi_mask.nii.gz}
    olsTD.nii.gz
   )

itk_add_test(NAME "${TestNamePrefix}itkRSHReconImageFilterTest"
  COMMAND ${itk-module}UnitTestDriver itkRSHReconImageFilterTest
    DATA{Baseline/template_dwi.nii.gz}
    DATA{Baseline/template_dwi.bvec}
    DATA{Baseline/template_dwi.bval}
    DATA{Baseline/template_dwi_mask.nii.gz}
    odf.nii.gz
   )

### itk_add_test(NAME "${TestNamePrefix}itkFiberTrackingManagerTest2"
###   COMMAND ${itk-module}UnitTestDriver itkFiberTrackingManagerTest2 )
