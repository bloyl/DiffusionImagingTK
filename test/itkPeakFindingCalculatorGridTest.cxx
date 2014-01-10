


#include <itkSymRealSphericalHarmonicRep.h>
#include <itkPeakFindingCalculatorGrid.h>
#include <itkPeakFindingCalculator.h>
#include <itkTimeProbe.h>

#include <itkTestingMacros.h>
namespace{


int test1()
{

  const unsigned int OdfOrder = 8;
  typedef double                                                          OdfPrecisionType;
  typedef itk::SymRealSphericalHarmonicRep<OdfPrecisionType,OdfOrder>     OdfPixelType;

  typedef itk::PeakFindingCalculatorGrid<OdfPixelType>                    PeakFinderType;
  typedef PeakFinderType::PeakDirectionContainerType                      PeakDirectionContainerType;
  typedef PeakFinderType::PeakValueContainerType                          PeakValueContainerType;

  //New lets test these 3 odfs.
  double ch[45] = {
         2.82094806e-01,   4.59387809e-01,   7.62639241e-03,
        -2.65089244e-01,  -1.70459377e-03,  -2.27291649e-03,
         3.38489264e-01,   8.24082829e-03,  -2.55495191e-01,
        -8.35784525e-03,   1.70923799e-01,   1.99015043e-03,
         1.02482527e-03,  -5.59562538e-03,  -3.80363222e-03,
         1.70554623e-01,   5.64062875e-03,  -1.26915440e-01,
        -5.26009873e-03,   1.14953570e-01,   4.45495034e-03,
        -7.90442452e-02,  -1.32022612e-03,  -3.45894368e-04,
         3.95073835e-03,   9.47112043e-04,  -6.37689978e-03,
        -3.68007110e-03,   4.89506535e-02,   2.26373016e-03,
        -3.73514183e-02,  -2.18456332e-03,   3.32550853e-02,
         1.17909722e-03,  -3.10148522e-02,  -6.25214074e-04,
         2.12846715e-02,   3.60342383e-04,   3.37901074e-05,
        -1.33388466e-03,  -1.39679192e-04,   2.47698487e-03,
         2.21211842e-04,  -3.72869126e-03,  -2.16668053e-03
     };
  double cc[45] = {
         2.82094806e-01,  -1.84499600e-03,  -4.37358674e-03,
        -2.57216394e-01,  -2.21597566e-03,  -5.69619704e-03,
         3.87160391e-01,  -5.81498304e-03,  -4.88071050e-03,
         3.55852279e-03,   1.55444667e-01,   1.25033641e-03,
         2.96612363e-03,   2.67166062e-03,  -6.83821831e-03,
        -1.60360534e-03,  -3.04495008e-03,  -1.28530160e-01,
         2.44779466e-03,   5.92144812e-03,  -4.25606035e-04,
        -6.48434013e-02,   5.08837402e-04,  -1.10527303e-03,
        -2.69610246e-05,   2.88396259e-03,  -2.25038873e-03,
        -7.32504530e-03,   1.14771903e-01,  -3.58927739e-03,
        -2.93563283e-03,   1.99614253e-04,   2.79063080e-02,
         2.29189231e-04,  -2.89570610e-03,  -1.19774416e-03,
         1.54047897e-02,  -7.50655425e-04,   2.04480690e-04,
        -1.01385545e-03,  -9.12652526e-04,  -5.92372544e-06,
         1.25703122e-03,   2.75207148e-03,  -3.13866278e-03
     };

  double cv[45] = {
         2.82094806e-01,  -4.56637770e-01,  -2.07055407e-03,
        -2.66447693e-01,  -3.83342733e-03,  -2.03245622e-03,
         3.29349637e-01,   6.67131925e-03,   2.56829888e-01,
         2.59662815e-03,   1.74307019e-01,   4.18567704e-03,
         2.51308928e-04,   3.66590638e-03,   3.40090320e-03,
        -1.61080286e-01,  -6.84451079e-03,  -1.24465309e-01,
        -5.27192699e-03,  -1.18184127e-01,  -1.89521688e-03,
        -8.28848481e-02,  -2.20265985e-03,   4.79890703e-04,
        -2.14970903e-03,   4.34971036e-04,  -1.83163770e-03,
        -3.47044342e-03,   4.55690846e-02,   3.02502280e-03,
         3.43799517e-02,   3.41011304e-03,   3.30946259e-02,
         2.10025162e-03,   3.35046090e-02,   6.53992000e-04,
         2.37989724e-02,   3.83429928e-04,  -1.72163273e-04,
         3.30030743e-04,  -7.50337844e-04,   3.65657092e-04,
        -9.34689073e-04,   3.83923936e-04,   2.20232457e-03
     };

  OdfPixelType odf_h(ch);
  OdfPixelType odf_c(cc);
  OdfPixelType odf_v(cv);
  OdfPixelType odf_i(0.0);
  odf_i[0] = 2.82094806e-01;

  unsigned int resol = 3;
  unsigned int rigid = 2;

  //compute peaks for each of these...
  PeakFinderType::Pointer finder = PeakFinderType::New();
  EXERCISE_BASIC_OBJECT_METHODS(finder, PeakFinderType);
  
  finder->Initialize(resol,rigid);

  PeakDirectionContainerType peaks;
  PeakValueContainerType peakVals;

  PeakDirectionContainerType::const_iterator peakIter;
  PeakValueContainerType::const_iterator valIter;

  //Look at the isotropic voxel... NO PEAKS...
  std::cout << "Testing isotropic" << std::endl;

  peaks = finder->ComputePeaksRaw(odf_i,peakVals);

  if (peakVals.size() != 0)
  {
    std::cout << "Isotropic Shouldn't have peaks"<< std::endl << "Raw Peaks: " << std::endl;
    valIter = peakVals.begin();
    for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
    {
      std::cout << *peakIter <<  " = " << *valIter << std::endl;
      ++valIter;
    }
    return EXIT_FAILURE;
  }
  std::cout << "\tComputePeaksRaw : passed" << std::endl;
  
  peaks = finder->ComputePeaks(odf_i,peakVals);
  if (peakVals.size() != 0)
  {
    std::cout << "Isotropic Shouldn't have peaks"<< std::endl << "Peaks: " << std::endl;
    valIter = peakVals.begin();
    for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
    {
      std::cout << *peakIter <<  " = " << *valIter << std::endl;
      ++valIter;
    }
    return EXIT_FAILURE;
  }
  std::cout << "\tComputePeaks    : passed" << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  //Look at the horizontal voxel...
  std::cout << "Testing Horizontal" << std::endl;

  peaks = finder->ComputePeaksRaw(odf_h,peakVals);

  if (peakVals.size() != 3)
  {
    std::cout << "Horizontal Should have 3 raw peaks"<< std::endl;
    return EXIT_FAILURE;
  }

  // valIter = peakVals.begin();
  // for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << "\tComputePeaksRaw : passed" << std::endl;
  
  // peaks = finder->ComputePeaks(odf_h,peakVals);
  // if (peakVals.size() != 0)
  // {
  //   std::cout << "Horizontal Shouldn't have peaks"<< std::endl << "Peaks: " << std::endl;
  //   valIter = peakVals.begin();
  //   for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  //   {
  //     std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //     ++valIter;
  //   }
  //   return EXIT_FAILURE;
  // }
  // std::cout << "\tComputePeaks    : passed" << std::endl;



  // std::cout << std::endl << "Raw Peaks: " << std::endl;
  // valIter = rawPeakVals.begin();
  // for (peakIter = rawPeaks.begin(); peakIter != rawPeaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;
  // std::cout << std::endl << "peaks: " << std::endl;
  // valIter = peakVals.begin();
  // for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;

  // std::cout << "Looking at Horizontal" << std::endl << odf_h << std::endl;
  // rawPeaks = finder->ComputePeaksRaw(odf_h,rawPeakVals);
  
  // std::cout << std::endl << "Raw Peaks: " << std::endl;
  // valIter = rawPeakVals.begin();
  // for (peakIter = rawPeaks.begin(); peakIter != rawPeaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;
  // std::cout << std::endl << "peaks: " << std::endl;
  
  // peaks = finder->ComputePeaks(odf_h,peakVals);
  // valIter = peakVals.begin();
  // for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;

  // return EXIT_FAILURE;
  
  // std::cout << "Looking at Vertical" << std::endl;
  // rawPeaks = finder->ComputePeaksRaw(odf_v,rawPeakVals);
  // peaks = finder->ComputePeaks(odf_v,peakVals);
  
  // std::cout << std::endl << "Raw Peaks: " << std::endl;
  // valIter = rawPeakVals.begin();
  // for (peakIter = rawPeaks.begin(); peakIter != rawPeaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;
  // std::cout << std::endl << "peaks: " << std::endl;
  // valIter = peakVals.begin();
  // for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;


  // std::cout << "Looking at crossing" << std::endl;
  // rawPeaks = finder->ComputePeaksRaw(odf_c,rawPeakVals);
  // peaks = finder->ComputePeaks(odf_c,peakVals);
  
  // std::cout << std::endl << "Raw Peaks: " << std::endl;
  // valIter = rawPeakVals.begin();
  // for (peakIter = rawPeaks.begin(); peakIter != rawPeaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;
  // std::cout << std::endl << "peaks: " << std::endl;
  // valIter = peakVals.begin();
  // for (peakIter = peaks.begin(); peakIter != peaks.end(); ++peakIter )
  // {
  //   std::cout << *peakIter <<  " = " << *valIter << std::endl;
  //   ++valIter;
  // }
  // std::cout << std::endl;

  // //Do timing tests...
  // itk::TimeProbe * oldCalcProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * initCalcProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * rawCalcProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * fullCalcProbe = new itk::TimeProbe(); 
  
  // int nReps =  200;
  // std::vector< double >            odPeaks; // each peak by 4 components: vx, vy, vz, mag

  // typedef itk::PeakFindingCalculator < OdfPixelType, double, double > PeakFindingCalculator;
  // PeakFindingCalculator::Pointer       pFinder = PeakFindingCalculator::New();

  // int nmPeaks;
  // oldCalcProbe->Start();
  // for (int j = 0; j < nReps; j++)
  // {
  //   for ( int i = 0; i < OdfPixelType::Dimension ; i ++ )
  //   {
  //     pFinder->SetCoefficient(  i,  odf_h.GetNthComponent( i )  );
  //   }
  //   nmPeaks   = pFinder->GetRawPeaks( odPeaks );
  // }  
  // oldCalcProbe->Stop();
  // std::cout << "Old Method ComputePeaks (" << nReps << " reps) elapsed time: " << oldCalcProbe->GetMean() << std::endl;

  // PeakFinderType::Pointer finder = PeakFinderType::New();; 
  // initCalcProbe->Start();
  // finder->Initialize(resol,rigid);
  // initCalcProbe->Stop();
  // std::cout << "Initialization elapsed time: " << initCalcProbe->GetMean() << std::endl;

  // rawPeaks;
  // rawCalcProbe->Start();
  // for (int i = 0; i < 200; i++)
  // {
  //   rawPeaks = finder->ComputePeaksRaw(odf_h);
  // }  
  // rawCalcProbe->Stop();
  // std::cout << "ComputeRawPeaks         (" << nReps << " reps) elapsed time: " << rawCalcProbe->GetMean() << std::endl;

  // peaks;
  // fullCalcProbe->Start();
  // for (int i = 0; i < 200; i++)
  // {
  //   peaks = finder->ComputePeaks(odf_h);
  // }  
  // fullCalcProbe->Stop();
  // std::cout << "ComputePeaks            (" << nReps << " reps) elapsed time: " << fullCalcProbe->GetMean() << std::endl;
  
  std::cout << "*************************************************************************************" << std::endl
            << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;
  
  return EXIT_FAILURE;
}

} // end empty namespace

using namespace itk;
int itkPeakFindingCalculatorGridTest( int, char **)
{

  if ( test1() == EXIT_FAILURE )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
