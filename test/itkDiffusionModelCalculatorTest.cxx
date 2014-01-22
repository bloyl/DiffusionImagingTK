
#include <itkSymRealSphericalHarmonicRep.h>
#include <itkDiffusionModelCalculator.h>
#include <itkPeakFindingCalculatorGrid.h>
#include <itkPeakFindingCalculator.h>

#include <itkTimeProbe.h>
#include <itkVariableLengthVector.h>
#include <vnl/vnl_random.h>

#include <itkTestingMacros.h>

namespace{

//GradientDirectionContainerType
template <class GradientDirectionContainerType>
typename GradientDirectionContainerType::Pointer generateGradientDirections( int resolution )
{
  typename GradientDirectionContainerType::Pointer gradCont = GradientDirectionContainerType::New();
  typename GradientDirectionContainerType::Element gradDir;
  
  if (resolution == -1)
  {
    //Use the BBL stuff
    const unsigned int numberOfGradientImages = 67; // The bbl set!!!

    //Set up Gradient Contatiner...
    // typedef typename GradientDirectionContainerType::Element       GradientDirectionType;
    
    // GradientDirectionType dir;
    double  gradientDirections[numberOfGradientImages][3] =
    {
      {0.000000,0.000000,0.000000},
      {1.000000,0.000000,0.000000},
      {0.000000,1.000000,0.000000},
      {-0.026007,0.649170,0.760199},
      {0.591136,-0.766176,0.252058},
      {-0.236071,-0.524158,0.818247},
      {-0.893021,-0.259006,0.368008},
      {0.796184,0.129030,0.591137},
      {0.233964,0.929855,0.283956},
      {0.935686,0.139953,0.323891},
      {0.505827,-0.844710,-0.174940},
      {0.346220,-0.847539,-0.402256},
      {0.456968,-0.630956,-0.626956},
      {-0.486997,-0.388997,0.781995},
      {-0.617845,0.672831,0.406898},
      {-0.576984,-0.104997,-0.809978},
      {-0.826695,-0.520808,0.212921},
      {0.893712,-0.039987,-0.446856},
      {0.290101,-0.541189,-0.789276},
      {0.115951,-0.962591,-0.244896},
      {-0.800182,0.403092,-0.444101},
      {0.513981,0.839970,0.173994},
      {-0.788548,0.152912,-0.595659},
      {0.949280,-0.233069,0.211062},
      {0.232964,0.782880,0.576911},
      {-0.020999,-0.187990,-0.981946},
      {0.216932,-0.955701,0.198938},
      {0.774003,-0.604002,0.190001},
      {-0.160928,0.355840,0.920587},
      {-0.147035,0.731173,-0.666158},
      {0.888141,0.417066,0.193031},
      {-0.561971,0.231988,-0.793959},
      {-0.380809,0.142928,0.913541},
      {-0.306000,-0.199000,-0.931001},
      {-0.332086,-0.130034,0.934243},
      {-0.963226,-0.265062,0.044010},
      {0,0,0},
      {-0.959501,0.205107,0.193101},
      {0.452965,-0.888932,0.067995},
      {-0.773133,0.628108,0.088015},
      {0.709082,0.408047,0.575066},
      {-0.692769,0.023992,0.720760},
      {0.681659,0.528735,-0.505747},
      {-0.141995,-0.724976,0.673978},
      {-0.740168,0.388088,0.549125},
      {-0.103006,0.822044,0.560030},
      {0.584037,-0.596038,0.551035},
      {-0.088008,-0.335031,0.938088},
      {-0.552263,-0.792377,0.259123},
      {0.838158,-0.458086,-0.296056},
      {0.362995,-0.560993,0.743990},
      {-0.184062,0.392133,-0.901306},
      {-0.720938,-0.692941,0.008999},
      {0.433101,0.682159,-0.589137},
      {0.502114,0.690157,0.521119},
      {-0.170944,-0.508833,-0.843722},
      {0.462968,0.422971,0.778946},
      {0,0,0},
      {0.385030,-0.809064,0.444035},
      {-0.713102,-0.247035,0.656094},
      {0.259923,0.884737,-0.386885},
      {0.001000,0.077002,-0.997030},
      {0.037002,-0.902057,0.430027},
      {0.570320,-0.303170,-0.763428},
      {-0.282105,0.145054,-0.948354},
      {0.721098,0.608082,0.332045},
      {0.266985,0.959945,-0.084995}
    };
    
    for (unsigned int g = 0; g<numberOfGradientImages;++g)
    {
      gradDir[0] = gradientDirections[g][0];
      gradDir[1] = gradientDirections[g][1];
      gradDir[2] = gradientDirections[g][2];
      gradCont->InsertElement(g,gradDir);
    }
  }
  else
  {
    typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double, double> MeshTraits;
    typedef itk::Mesh<double,3,MeshTraits> TriangleMeshType;

    // declare triangle mesh source
    typedef itk::RegularSphereMeshSource<TriangleMeshType>  SphereMeshSourceType;
    typedef SphereMeshSourceType::PointType PointType;
    typedef SphereMeshSourceType::VectorType VectorType;

    SphereMeshSourceType::Pointer  mySphereMeshSource = SphereMeshSourceType::New();
    PointType center; center.Fill(0);
    PointType::ValueType scaleInit[3] = {1,1,1};
    VectorType scale = scaleInit;

    mySphereMeshSource->SetCenter(center);
    mySphereMeshSource->SetResolution(resolution); 
    mySphereMeshSource->SetScale(scale);
    mySphereMeshSource->Update();
    
    TriangleMeshType::Pointer sphere = mySphereMeshSource->GetOutput();

    unsigned int numPoints = sphere->GetPoints()->Size();
    PointType  point(0);

    gradCont->Reserve(numPoints+resolution);
    for (unsigned int pointIndex = 0; pointIndex < numPoints; pointIndex++)
    {
      sphere->GetPoint(pointIndex,&point);
      gradDir[0] = point[0]; gradDir[1] = point[1]; gradDir[2] = point[2];
      gradCont->InsertElement(pointIndex,gradDir);
    }

    //Add 0,0,0 vectors for B0images
    for (int pointIndex = 0; pointIndex < resolution; pointIndex++)
    {
      gradDir[0] = 0.0; gradDir[1] = 0.0; gradDir[2] = 0.0;
      gradCont->InsertElement(numPoints+pointIndex,gradDir);
    }
  }

  return gradCont;
}

template <typename DWIPixelType, typename GradientDirectionContainerType, typename DtType>
DWIPixelType
generateTestData(typename GradientDirectionContainerType::Pointer gradContainer,
                DtType dt, double s0, double noiseAmp)

{
  const unsigned int numberOfGradientImages = gradContainer->Size();
  typename GradientDirectionContainerType::Element dir;

  DWIPixelType dwi;
  dwi.SetSize(numberOfGradientImages);
  dwi.Fill(0);

  //Convert the dt to vnl mats
  vnl_matrix_fixed<double,3,3> T;
  T[0][0] = dt[0];
  T[0][1] = T[1][0] = dt[1];
  T[0][2] = T[2][0] = dt[2];
  T[1][1] = dt[3];
  T[1][2] = T[2][1] = dt[4];
  T[2][2] = dt[5];

  vnl_random randGen;
  for (unsigned int i=0;i<numberOfGradientImages;++i)
  {
    dir = gradContainer->ElementAt(i);
    dwi[i] = s0 * exp(-dot_product(dir, T * dir));
    if (noiseAmp > 0)
    {
      dwi[i] += abs(noiseAmp * randGen.normal());
    }
  }
  return dwi;
}

//TODO add test for SPD...itkDiffusionModelCalculatorTest

int diffusionTensorTest()
{
  typedef itk::VariableLengthVector<unsigned int>                 DWIPixelType;

  typedef itk::DiffusionModelCalculator<DWIPixelType, double>     DTCalculatorType;
  typedef DTCalculatorType::GradientDirectionContainerType        GradientDirectionContainerType;
  typedef DTCalculatorType::DtType                                DtType;
  typedef DTCalculatorType::RshType                               RshType;

  GradientDirectionContainerType::Pointer gradDirs = generateGradientDirections<GradientDirectionContainerType>(3);
  double bValue = 1000.0;
  double bValueSqrt = vcl_sqrt(bValue);

  std::cout << "Using " << gradDirs->Size() << " gradient Directions" << std::endl;
  
  //Gradient needs to be scaled by the sqrt of the bvalue...
  for (unsigned int i=0;i<gradDirs->Size();++i)
  {
    gradDirs->SetElement(i, bValueSqrt * gradDirs->ElementAt(i) );
  }

  DtType inpDt1;
  DtType inpDt2;
  DtType inpDt3;
  
  inpDt1[0] = 0.2E-3;
  inpDt1[1] = 0;
  inpDt1[2] = 0;
  inpDt1[3] = 0.5E-3;
  inpDt1[4] = 0;
  inpDt1[5] = 1.7E-3;
  inpDt2 = inpDt1;

  inpDt3[0] = -0.2E-3;
  inpDt3[1] = 0;
  inpDt3[2] = 0;
  inpDt3[3] = 0.5E-3;
  inpDt3[4] = 0;
  inpDt3[5] = 1.7E-3;

  //Test No Noise!
  DWIPixelType sig1 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt1, 5000., 0.);
  DWIPixelType sig2 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt2, 5000., 0.);
  DWIPixelType sig3 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt3, 5000., 0.);

  DWIPixelType res1;
  res1.SetSize(sig1.Size());

  DTCalculatorType::Pointer dtCalc = DTCalculatorType::New();
  EXERCISE_BASIC_OBJECT_METHODS(dtCalc, DTCalculatorType);

  TRY_EXPECT_EXCEPTION(dtCalc->InitializeTensorFitting());

  GradientDirectionContainerType::Pointer gradDirsSmall = GradientDirectionContainerType::New();
  gradDirsSmall->Reserve(4);
  for (unsigned int i=0;i<4;++i)
  {
    gradDirsSmall->InsertElement(i, gradDirs->ElementAt(i));
  }

  dtCalc->SetGradientDirectionContainer(gradDirsSmall);
  TRY_EXPECT_EXCEPTION(dtCalc->InitializeTensorFitting());

  dtCalc->SetGradientDirectionContainer(gradDirs);
  TEST_SET_GET(gradDirs, dtCalc->GetGradientDirectionContainer());

  dtCalc->InitializeTensorFitting();

  //Test ordinary least squares (OLS)
  DtType dt1    = dtCalc->ComputeTensorOLS(sig1, false);
  DtType dt2    = dtCalc->ComputeTensorOLS(sig2, false);
  DtType dt3    = dtCalc->ComputeTensorOLS(sig3, false);
  DtType dt3spd = dtCalc->ComputeTensorOLS(sig3, true);

  std::cout << "dt ols fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares (WLS)
  dt1    = dtCalc->ComputeTensorWLS(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS(sig3, true);

  std::cout << "dt WLS fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares with residual boot straping(WLS)
  dt1    = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS_residualBoot(sig3, true);

  std::cout << "dt WLS residual boot strap fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  // Test with noise
  sig1 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt1, 5000., 500.);
  sig2 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt2, 5000., 500.);
  sig3 = generateTestData<DWIPixelType, GradientDirectionContainerType, DtType>(gradDirs, inpDt3, 5000., 500.);

  //Test ordinary least squares (OLS)
  dt1    = dtCalc->ComputeTensorOLS(sig1, false);
  dt2    = dtCalc->ComputeTensorOLS(sig2, false);
  dt3    = dtCalc->ComputeTensorOLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorOLS(sig3, true);

  std::cout << "dt ols fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares (WLS)
  dt1    = dtCalc->ComputeTensorWLS(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS(sig3, true);

  std::cout << "dt WLS fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;

  //Test weighted least squares with residual boot straping(WLS)
  dt1    = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  dt2    = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  dt3    = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);
  dt3spd = dtCalc->ComputeTensorWLS_residualBoot(sig3, true);

  std::cout << "dt WLS residual boot strap fit" << std::endl;
  std::cout << "sig1      : " << dt1 << std::endl;
  std::cout << "sig2      : " << dt2 << std::endl;
  std::cout << "sig3      : " << dt3 << std::endl;
  std::cout << "sig3 spd  : " << dt3spd << std::endl;
  std::cout << "!!!! NO CHECKS Made !!!" << std::endl;


  // dtCalc->InitializeRSHFitting();

  // RshType adc1  = dtCalc->ComputeRSH_ADC(sig1);
  // RshType odf1  = dtCalc->ComputeRSH_ODF(sig1);
  // RshType csaodf = dtCalc->ComputeRSH_CSAODF(sig1);
  // RshType csaodf_noReg = dtCalc->ComputeRSH_CSAODF(sig1,0,0);

  // std::cout << "adc1         " << adc1 << std::endl;
  // std::cout << "odf1         " << odf1 << std::endl;
  // std::cout << "csaodf       " << csaodf << std::endl;
  // std::cout << "csaodf_noReg " << csaodf_noReg << std::endl;

  // //Test Some Noise!
  // generateTestData<DWIPixelType, GradientDirectionContainerType> (gradDirs, bValue, inpDt1, inpDt2, inpDt3,5000,500, sig1, sig2, sig3);

  // dt1 = dtCalc->ComputeTensorOLS(sig1, false);
  // dt2 = dtCalc->ComputeTensorOLS(sig2, false);
  // dt3 = dtCalc->ComputeTensorOLS(sig3, false);

  // wlsDt1 = dtCalc->ComputeTensorWLS(sig1, false);
  // wlsDt2 = dtCalc->ComputeTensorWLS(sig2, false);
  // wlsDt3 = dtCalc->ComputeTensorWLS(sig3, false);

  // resBootDt1 = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  // resBootDt2 = dtCalc->ComputeTensorWLS_residualBoot(sig2, false);
  // resBootDt3 = dtCalc->ComputeTensorWLS_residualBoot(sig3, false);

  // std::cout << "sig1 " << std::endl;
  // std::cout << "ols         " << dt1 << std::endl;
  // std::cout << "wls         " << wlsDt1 << std::endl;
  // std::cout << "resBootDt   " << resBootDt1 << std::endl;

  // std::cout << "sig2 " << std::endl;
  // std::cout << "ols         " << dt2 << std::endl;
  // std::cout << "wls         " << wlsDt2 << std::endl;
  // std::cout << "resBootDt   " << resBootDt2 << std::endl;

  // std::cout << "sig3 " << std::endl;
  // std::cout << "ols         " << dt3 << std::endl;
  // std::cout << "wls         " << wlsDt3 << std::endl;
  // std::cout << "resBootDt   " << resBootDt3 << std::endl;

  // //Do some timing tests...
  // itk::TimeProbe * olsProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * wlsProbe = new itk::TimeProbe(); 
  // itk::TimeProbe * resBootProbe = new itk::TimeProbe(); 
  // olsProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorOLS(sig1, false);
  // olsProbe->Stop();

  // wlsProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorWLS(sig1, false);
  // wlsProbe->Stop();

  // resBootProbe->Start();
  // for (int i = 0; i < 200; i++)
  //   dt1 = dtCalc->ComputeTensorWLS_residualBoot(sig1, false);
  // resBootProbe->Stop();

  // std::cout << "200 ols elapsed time            : " << olsProbe->GetMean() << std::endl;
  // std::cout << "200 wls elapsed time            : " << wlsProbe->GetMean() << std::endl;
  // std::cout << "200 res Bootstrap elapsed time  : " << resBootProbe->GetMean() << std::endl;

  // std::cout << "*************************************************************************************" << std::endl
  //           << "***  NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
  //           << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}

} // end empty namespace

using namespace itk;
int itkDiffusionModelCalculatorTest( int , char ** )
{

  if ( diffusionTensorTest() == EXIT_FAILURE )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
