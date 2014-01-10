/**
 * @file  itkSymRealSphericalHarmonicRepTest.cxx
 * @brief Test itkSymRealSphericalHarmonicRep module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkSymRealSphericalHarmonicRep.h>
#include <limits.h>
#include <cmath>
#include <stdio.h>

#include <itkMesh.h>
#include <itkRegularSphereMeshSource.h>
#include <itkDefaultDynamicMeshTraits.h>

#include <itkTestingMacros.h>

namespace
{ //Empty Namespace

//We currently only suporrt Orders up to 20 so only test that HIGH!
const int maxOrder = 20;
typedef itk::SymRealSphericalHarmonicRep<double,maxOrder>   PixelType;

const int maxJ = PixelType::GetJ(maxOrder,maxOrder);

const double th1 = 1.324;
const double th2 = 0.898;
const double ph1 = 5.231;
const double ph2 = 1.321;
double percision  = 1e-6;
bool passed = true;

bool areEqual(double x, double y)
{
  //Compare binanry significant and eponent
  // x = xSig * 2^xExp
  double xSig,ySig;
  int xExp,yExp;

  xSig = frexp(x , &xExp);
  ySig = frexp(y , &yExp);

  if (xExp != yExp)
  {
    std::cerr << "areEqual exponents differ : " << xExp << " : " << yExp << std::endl;
    return false;
  }

  if ( vcl_abs(xSig - ySig) > percision )
  {
    std::cerr << xSig << " : " << ySig << std::endl;
    std::cerr << vcl_abs(xSig - ySig) << std::endl;
    return false;
  }
  return true;
}

//GradientDirectionContainerType
template <class GradientDirectionContainerType>
typename GradientDirectionContainerType::Pointer generateGradientDirections(unsigned int resolution)
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

  typename GradientDirectionContainerType::Pointer gradCont = GradientDirectionContainerType::New();
  typename GradientDirectionContainerType::Element gradDir;

  gradCont->Reserve(numPoints);
  for (unsigned int pointIndex = 0; pointIndex < numPoints; pointIndex++)
  {
    sphere->GetPoint(pointIndex,&point);
    gradDir[0] = point[0]; gradDir[1] = point[1]; gradDir[2] = point[2];
    gradDir.normalize();
    gradCont->InsertElement(pointIndex,gradDir);
  }

  return gradCont;
}

//------------------------------------------------------------------------------------------------------------

int testIndexing()
{
  typedef PixelType::LmVector                      LmVectorType;

  //test GetLM and GetJ
  int c = 1;
  for( int l=0;l<=maxOrder;l=l+2 )
  {
    for (int m=-l;m<=l;m++)
    {
      int j = PixelType::GetJ(l,m);
      if (j!=c)
      {
        std::cerr << "GetJ failed for " << l << " , " << m <<std::endl;
        std::cerr << "       Expected" << c << " but recieved " << j <<std::endl;
        passed = false;
      }
      LmVectorType vec = PixelType::GetLM(j);
      if ( vec(0) != l || vec(1) != m )
      {
        std::cerr << "GetLM failed for j = " << j << std::endl;
        std::cerr << "       Expected" << l << " , " << m << " but recieved " << vec <<std::endl;
        passed = false;
      }
      c++;
    }
  }
  //~ std::cerr << "last J = " << PixelType::GetJ(20,20)<<std::endl;

  if (!passed)
    return EXIT_FAILURE;

  std::cerr << "Passed GetLM / GetJ Testing" << std::endl;
  return EXIT_SUCCESS;
}


int testBasisFunctions()
{
  if (231 != maxJ)
    return EXIT_FAILURE;

  const double expResults1[231] = {0.28209479177387814,-0.26126056996555863,0.12829079677164226,-0.25892214317241713,
   0.2247880371192935,-0.44226847719203827,-0.2670659586656951,-0.39426660282468,0.1317337542965475,
   -0.20286428986633404,0.14114042500034443,-0.35545391149002886,0.223002219260276,-0.00589995716023252,
   0.48465065702649496,0.5677638731506147,0.2584882063484731,0.07396087053703475,0.48084600806007466,
   -0.009531982668877918,0.2276174825965264,0.01231287085239147,0.39882585824128197,-0.01613598048929191,
   0.007195564696319936,-0.13421847051984379,-0.4229919636698713,-0.016996280176664557,
   -0.30440068056043845,0.2696061245072423,-0.04635127841660211,-0.2643651274210245,0.06537157888489939,
   -0.4366902190802431,-0.10701669634532471,-0.19764234604319275,-0.16344262542133395,
   -0.3463041475826414,-0.18116056063495897,-0.0065348004786793205,-0.11863128799984893,
   0.43260899966516203,0.0013875474502164549,0.5070648010222,-0.48177077886920183,-0.2573648960487785,
   -0.6349712734191281,-0.031143054475546654,-0.23530221855449435,-0.26621538988912913,
   0.19401450919115626,-0.16887516734377164,0.2852936091587839,0.19426985419169226,0.12042594141299212,
   0.27587268361433137,0.2110074274201813,0.3288649052132344,0.004269243349717421,0.3064615992900459,
   -0.3174867410104502,0.00796928365446623,-0.44254733771158017,-0.04928968484376007,
   -0.028522897264176582,0.5020224166701597,0.5525945989942991,0.37355445663928205,-0.07229879616427291,
   0.4716349303990548,0.19980957249903517,0.13234760456648478,0.44327715030567766,-0.07754416818699836,
   0.22151763397645138,-0.06742521065385616,-0.23296332732739655,-0.014551226583123687,
   -0.3226346352574463,-0.02549630794733262,-0.39436619169995896,-0.001008976797748955,
   -0.4019930783617018,0.1268938356450434,-0.01326971122819245,0.24891427039324066,0.3162358677751695,
   0.021185831909422043,0.14102784385113826,-0.5720803885852106,-0.03311403612080456,
   -0.30188852771294816,0.31869804660523476,0.2523107602237113,-0.2327801957439004,0.20949343524183375,
   -0.17726539557007787,-0.2653206498482737,-0.003112507837851088,-0.47716023233248295,
   -0.05275092010147502,-0.2153234123021843,-0.16298427995700066,0.21478306732207295,
   -0.09464131470256222,0.2924850926635978,-0.16582822694792712,0.3635901893794796,-0.002438959482361345,
   0.39075228369359816,0.08632198580477797,0.014284017320231735,-0.005853884700744638,
   -0.4199193506799875,-0.00796275812465226,-0.4086431453989201,0.356491489980281,-0.015119633168549718,
   0.6452647332229915,-0.44755326745171875,-0.22500045471758875,-0.7447217436906549,-0.19033404566223933,
   -0.16241574219194913,-0.5219237568842503,0.03644690208488729,-0.2349850257961194,-0.12935740840179,
   0.24028395893913762,-0.11641698560327764,0.37809771906069234,0.16439086481464257,0.15475459095768307,
   0.350415437241157,-0.14471851808638975,0.1810693689440262,-0.19256404821785997,0.3172653771868516,
   -0.24498315464899123,0.005243751444314719,-0.28083666881481006,-0.26901039586924574,
   -0.011318534114636406,-0.21895257664647214,0.3802941236358348,-0.005810732271929013,
   0.4583676808398499,-0.05581664878270617,0.03127609674293864,-0.3288415215375242,-0.2821724452028133,
   -0.05582154729737142,0.47326216405778126,0.5047989414583224,0.437883242607259,-0.18071875901856832,
   0.2978706640117816,0.3167819823604304,-0.019751209740254565,0.5030787561808104,0.13791837565189444,
   0.1768684000462091,0.3647289782995581,-0.14624666014699744,0.19717371872083359,-0.1817217133182147,
   -0.23097553691292091,-0.056355120673604465,-0.450459728013501,0.03996183951188233,
   -0.22411726648813624,0.04668381338591086,-0.39269286407258197,0.06764840905400585,
   -0.006740852708925628,0.1022689166291643,0.377970033134846,0.005439925468664964,0.37083672573358395,
   -0.2314625816089435,0.01638361861834898,-0.3450039340522714,-0.21121524996867044,-0.0301468167334597,
   -0.03999007593436368,0.46963298787575125,0.022327267198122693,0.3801207916949088,-0.6286674578345203,
   -0.045442812001185294,-0.2853539187841342,0.3238583546873271,0.4922296638450851,-0.12827944234390604,
   0.2540929317748696,0.1516556293402713,-0.26202438128834693,0.15813918891286163,-0.30039162653143237,
   -0.24555161249761773,-0.06782596005621196,-0.47605613178882705,0.01545748587804534,-0.221615168852733,
   -0.0572139494330796,0.23776486913112466,-0.054998453468957716,0.4400650899310798,0.07416858050603845,
   0.21354432658798272,0.11035768154651847,0.3741672140139402,0.12555444229577004,0.006585303344756739,
   0.09980694185939816,-0.38908014530419405,0.001712726646699089,-0.4168052625035737,0.02446435072718901,
   -0.02138443218994445,0.13230301763430743,0.37605029038914745,0.01800086209574355,0.3201828270769258,
   -0.38845420488823995,0.011367536879179768,-0.5344547899447447,0.18417035196550002,-0.0443113054296403,
   0.7083401429422193,-0.3968846930273987};

  const double expResults2[231] = {0.28209479177387814,-0.29329185613168046,0.13163258988127263,0.052053061163180134,
   -0.5159533571757736,0.16007097856985705,0.1266651551584488,-0.35945288290398447,-0.43647885023047944,
   -0.022696611741745316,-0.3567028152343802,0.08896271839087028,0.2382186727001045,0.38631913720970934,
   -0.19691624873045233,-0.011249149440558728,0.4093666203104718,0.3341209039399373,-0.237913850587078,
   0.25060893306578813,-0.11747707772160564,0.1143856909836888,0.4604687387548323,-0.13677576214791934,
   0.25569602543315934,-0.5194312117170425,-0.13648471420478103,0.15592820114859382,
   -0.042305785480841906,-0.32005250390551016,-0.04229428087058343,0.5619813171411655,
   0.06834752343210033,0.3048454418865678,0.3369065998531418,0.0712765877765214,0.3058947832228189,
   -0.2793791019723905,-0.183874758189409,-0.32763022274435066,-0.10625446206914599,-0.1873671561210002,
   0.5862550915401387,-0.057505381052407996,-0.09283929370665928,0.052552190593683455,
   0.18242212792045578,-0.20610273991122174,-0.6554368850160447,-0.032951974792404394,
   -0.11801116789950078,-0.2900984455811412,0.15400736689323644,-0.3638831237535335,0.08627494230866085,
   -0.2529569368323972,-0.3381673654256099,0.1985978351820756,-0.16551819704844753,0.4509929948367749,
   0.03934546620220314,0.4567582803325459,-0.11776551462249649,-0.45228879660041094,0.14672608774396023,
   0.039424677129634295,-0.041499568188732704,-0.06281372678105904,0.3142077189002391,0.4939210992997543,
   -0.2669464286508061,-0.2360141274730389,0.025357647663271005,-0.5011581085002027,
   -0.011767888044058324,-0.3418485174527467,-0.19066228825845893,-0.10864512070141834,
   -0.19225467826791115,0.42585058014298777,0.10405846060791142,0.3673989849569223,0.01829460017016635,
   0.1670884185871886,-0.35149078660451183,-0.04240579957498353,-0.5858091892580257,0.39727148992609984,
   0.23571877269697955,-0.1511606945177426,-0.006019145431062305,0.024936400772157124,
   -0.01183943639221538,-0.2939651297098082,-0.21100191278060623,0.5612090758069225,0.4041698216543932,
   0.008687216637243493,0.5044446428368407,0.03029707020337359,0.18505059900600246,0.282709549877772,
   -0.031099766791204003,0.43849913543745245,-0.03846885154350947,0.3394525474077542,0.1507843393377915,
   -0.23932129121237977,0.033424227890836926,-0.43950606599396946,-0.061696720899233774,
   -0.4199577649712567,0.09063600831800572,0.019063942383096857,0.3250825839581352,0.42101930225850115,
   -0.5077742925787732,-0.04263704285638409,0.11152899878674122,-0.009259822710945017,
   -0.011019589138870868,0.043080012879933856,0.2041436896174781,-0.047428882276277826,
   -0.6654653373792453,-0.26010353303299966,0.24674129783487067,-0.2143865186334627,0.2377183845956529,
   0.2288902758331066,-0.02766691042071162,0.4392778626043064,-0.07809949485982122,0.34941162213892996,
   0.006701166450282706,0.12549375296141713,0.04006404259009357,-0.49189128014142486,
   -0.0036573203409188777,-0.3755273717802197,0.12141507691125601,-0.14645725997045136,
   0.3835002455995976,0.041125822702076216,0.5216687664844939,-0.17243574287134664,0.18510543312833733,
   -0.6259369203936915,-0.09651986321400177,0.4467861119914072,-0.07580622362817883,-0.06225227649622372,
   0.012674980561662962,0.0022652726202464966,-0.045135091354123524,-0.10235573700797508,
   0.20031375908221288,0.5580367795392146,-0.07708384130956969,-0.5513787846497726,-0.01842477784798096,
   -0.37981445393574326,-0.40732542780768227,-0.0010106051079454255,-0.5017036032523858,
   -0.023452026597667736,-0.30680516088980786,-0.2512503975564898,-0.10660157477768581,
   -0.43983691292945276,-0.0170463889180076,-0.3574813189647205,0.0668158363966157,0.24005141497059654,
   0.11456919766675598,0.39059902241440303,0.10229025187456428,0.3250763393255871,-0.09014351248104036,
   -0.0022177549328444316,-0.32762070666622456,-0.2849369749656992,-0.04433906979510322,
   -0.07997261748160649,0.7261396031947156,-0.2072200271375747,-0.28946109071829523,0.11773188279537847,
   0.02269201982594112,-0.010321047572507474,0.0018544837059451065,0.033387803068027994,
   0.023535616812975637,-0.2388343155445499,-0.3273636387058059,0.40167340036288784,0.6584975321700367,
   -0.03724369729746505,0.21484394426601727,0.22349459414813946,-0.2964799472337818,0.18651067466287113,
   -0.22597166641750652,-0.12806751676957653,0.03310738092020865,-0.33286169976355334,
   0.17015479509479617,-0.3042730462818312,0.18194901075069125,-0.11788037755310135,0.1200024802430867,
   0.4620495319474324,-0.09930298299044596,0.32701504510443863,-0.26452613516047235,0.11097762178921863,
   -0.45891241634825364,-0.02301051003320508,-0.4958908023919907,0.15001459487251703,-0.2224193903823256,
   0.5378378230947594,0.03116121448876437,0.35084037221855396,-0.24452487988152763,-0.580433521462941,
   0.3765410583734806,0.12007581819055889,-0.10723310678102649,0.001126973207440763,0.00636430374593853};

  passed = true;
  for (int i = 0;i < maxJ; i++)
  {
    double res1 = PixelType::Y(i+1,th1,ph1);
    if (vcl_fabs(res1 - expResults1[i]) >= percision){
      printf("Failed - %d : (%f,%f) || %f :: %f :: %f \n",i+1,th1,ph1,res1,expResults1[i],res1-expResults1[i]);
      passed = false;
    }
  }

  if (!passed)
  {
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  passed = true;
  for (int i = 0;i < maxJ; i++)
  {
    double res1 = PixelType::Y(i+1,th2,ph2);
    if (vcl_fabs(res1 - expResults2[i]) >= percision){
      printf("Failed - %d : (%f,%f) || %f :: %f :: %f \n",i+1,th2,ph2,res1,expResults2[i],res1-expResults2[i]);
      passed = false;
    }
  }

  if (!passed)
  {
    std::cout << "[FAILED]" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Passed Basis computation!" << std::endl;
  return EXIT_SUCCESS;

}

int testNormalize()
{
  std::cout << "Testing Normalization" << std::endl;

  typedef itk::SymRealSphericalHarmonicRep<double,8>    PixelType2;

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

  double newVals[45];
  for (int i=0;i<45;++i)
  {
    newVals[i]=2.1*ch[i];
  }
  
  PixelType2 pix1(newVals);
  PixelType2 pixN(newVals);
  pixN.Normalize();

  std::cout << pixN << std::endl;
  std::cout << "Name of Class = " << pixN.GetNameOfClass() << std::endl;
  
  //GetNameOfClass’ is not a member of ‘itk::FixedArray
  //std::cout << "Name of Superclass = " << PixelType2::Superclass::GetNameOfClass() << std::endl;

  for (int i=0;i<45;++i)
  {
    if (  vcl_abs( (pix1[i]/pixN[i])-2.1) > percision )
      {
        std::cout << i << " " << pix1[i] << " : " << pixN[i] << " : " << vcl_abs( (pix1[i]/pixN[i])) << std::endl;
        std::cout << "----- Expected Ratio 2.1" << std::endl;
        std::cout << "----- Normalization Failed" << std::endl;
        return EXIT_FAILURE;
      }
  }
  std::cout << "----- Normalization Passed" << std::endl;
  return EXIT_SUCCESS;
}

int testEvaluation()
{
  typedef itk::SymRealSphericalHarmonicRep<double,8>    PixelType2;
  typedef PixelType2::GradientDirectionType             GradType;
  typedef PixelType2::GradientDirectionContainerType    GradientDirectionContainerType;

  typedef vnl_vector_fixed< double, PixelType2::Length>
                                                        VnlVectorPixelType;

  GradType gradX;
  GradType gradY;
  GradType gradZ;
  gradX.fill(0);gradX[0]=1;
  gradY.fill(0);gradY[1]=1;
  gradZ.fill(0);gradZ[2]=1;

  GradientDirectionContainerType::Pointer grads = GradientDirectionContainerType::New();
  grads->InsertElement(0,gradX);
  grads->InsertElement(1,gradY);
  grads->InsertElement(2,gradZ);

  PixelType2::RshBasisMatrixType basisMat
      = PixelType2::ComputeRshBasis( grads );

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
  PixelType2 odf_h(ch);
  VnlVectorPixelType vnl_h(ch);

  PixelType2 odf_c(cc);
  VnlVectorPixelType vnl_c(cc);

  PixelType2 odf_v(cv);
  VnlVectorPixelType vnl_v(cv);

  PixelType2 odf_i(0.0);
  VnlVectorPixelType vnl_i(0.0);

  odf_i[0] = 2.82094806e-01;
  vnl_i[0] = 2.82094806e-01;

  typedef vnl_vector_fixed< double, 3>                      ValuesType;
  ValuesType evalVals;
  ValuesType basisVals;

  std::cout << "Testing Evaluate" << std::endl;

  //First test odf_i
  evalVals[0] = odf_i.Evaluate(gradX);
  evalVals[1] = odf_i.Evaluate(gradY);
  evalVals[2] = odf_i.Evaluate(gradZ);
  basisVals   = basisMat * vnl_i;

  for (unsigned int i =0; i<3; i++)
  {
    if ( vcl_abs(evalVals[i] - basisVals[i]) > 0.0001 )
    {
      std::cout << "Isotropic Pixels" << std::endl;
      std::cout << "evaluate returned different results then using a basis" << std::endl;
      if (i == 0)
        std::cout << "For X gradient" << std::endl;
      else if (i == 1)
        std::cout << "For Y gradient" << std::endl;
      if (i == 2)
        std::cout << "For Z gradient" << std::endl;
      std::cout << "evaluate : " << evalVals[i] << std::endl;
      std::cout << "basis : " << basisVals[i] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "IsoPixelPassed" << std::endl;

  evalVals[0] = odf_h.Evaluate(gradX);
  evalVals[1] = odf_h.Evaluate(gradY);
  evalVals[2] = odf_h.Evaluate(gradZ);
  basisVals   = basisMat * vnl_h;

  for (unsigned int i =0; i<3; i++)
  {
    if ( vcl_abs(evalVals[i] - basisVals[i]) > 0.0001 )
    {
      std::cout << "Horizontal Pixels" << std::endl;
      std::cout << "evaluate returned different results then using a basis" << std::endl;
      if (i == 0)
        std::cout << "For X gradient" << std::endl;
      else if (i == 1)
        std::cout << "For Y gradient" << std::endl;
      if (i == 2)
        std::cout << "For Z gradient" << std::endl;
      std::cout << "evaluate : " << evalVals[i] << std::endl;
      std::cout << "basis : " << basisVals[i] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "HorizontalPassed" << std::endl;

  evalVals[0] = odf_v.Evaluate(gradX);
  evalVals[1] = odf_v.Evaluate(gradY);
  evalVals[2] = odf_v.Evaluate(gradZ);
  basisVals   = basisMat * vnl_v;

  for (unsigned int i =0; i<3; i++)
  {
    if ( vcl_abs(evalVals[i] - basisVals[i]) > 0.0001 )
    {
      std::cout << "Vertical Pixels" << std::endl;
      std::cout << "evaluate returned different results then using a basis" << std::endl;
      if (i == 0)
        std::cout << "For X gradient" << std::endl;
      else if (i == 1)
        std::cout << "For Y gradient" << std::endl;
      if (i == 2)
        std::cout << "For Z gradient" << std::endl;
      std::cout << "evaluate : " << evalVals[i] << std::endl;
      std::cout << "basis : " << basisVals[i] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "VerticalPixelPassed" << std::endl;

  evalVals[0] = odf_c.Evaluate(gradX);
  evalVals[1] = odf_c.Evaluate(gradY);
  evalVals[2] = odf_c.Evaluate(gradZ);
  basisVals   = basisMat * vnl_c;
  for (unsigned int i =0; i<3; i++)
  {
    if ( vcl_abs(evalVals[i] - basisVals[i]) > 0.0001 )
    {
      std::cout << "Vertical Pixels" << std::endl;
      std::cout << "evaluate returned different results then using a basis" << std::endl;
      if (i == 0)
        std::cout << "For X gradient" << std::endl;
      else if (i == 1)
        std::cout << "For Y gradient" << std::endl;
      if (i == 2)
        std::cout << "For Z gradient" << std::endl;
      std::cout << "evaluate : " << evalVals[i] << std::endl;
      std::cout << "basis : " << basisVals[i] << std::endl;
      return EXIT_FAILURE;
    }
  }
  std::cout << "CrossingPixelPassed" << std::endl;

  //Test Component Access
  TRY_EXPECT_EXCEPTION( odf_c.GetLthMthComponent(22,0) );
  TRY_EXPECT_EXCEPTION( odf_c.GetLthMthComponent(11,0) );

  TRY_EXPECT_EXCEPTION( odf_c.GetLthMthComponent(4,-7) );

  unsigned int indexI = 0;
  for( int l=0;l<=PixelType2::MaxOrder;l=l+2 )
  {
    for (int m=-l;m<=l;m++)
    {
      double getNValue  = odf_c.GetNthComponent(indexI);
      double getLMValue = odf_c.GetLthMthComponent(l,m);
      bool pass = true;
      if ( getNValue != cc[indexI] )
      {
        std::cout << "GetNthComponent Failed\n";
        pass = false;
      }
      if ( getLMValue != cc[indexI] )
      {
        std::cout << "GetLthMthComponent Failed\n";
        pass = false;
      }

      if (! pass)
      {
        std::cout << "indexI " << indexI << std::endl;
        std::cout << "l " << l << std::endl;
        std::cout << "m " << m << std::endl;
        std::cout << "expected           : " << cc[indexI] << std::endl;
        std::cout << "GetNthComponent    : " << getNValue << std::endl;
        std::cout << "GetLthMthComponent : " << getLMValue << std::endl;
        return EXIT_FAILURE;
      }
      indexI++;
    }
  }

  std::cout << "Passed GetNthComponent and GetLthMthComponent" << std::endl;


  ///TODO Should test SetNthComponent and SetLthMthComponent

  return EXIT_SUCCESS;
}

int testCurvatureComputation()
{


  typedef itk::SymRealSphericalHarmonicRep<double,8>    PixelType2;
  typedef PixelType2::GradientDirectionType             GradType;
  typedef PixelType2::GradientDirectionContainerType    GradientDirectionContainerType;
  typedef vnl_vector_fixed< double, PixelType2::Length>
                                                        VnlVectorPixelType;

  GradientDirectionContainerType::Pointer grads = generateGradientDirections<GradientDirectionContainerType>(5);

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

  PixelType2 odf_h(ch);
  PixelType2 odf_c(cc);
  PixelType2 odf_v(cv);
  PixelType2 odf_i(0.0);
  odf_i[0] = 2.82094806e-01;

  //first investigate the theta dependence...
  unsigned nPnts = grads->Size();

  typedef vnl_vector< double>                      ValuesType;
  ValuesType k1CurvaturesDir(nPnts);
  ValuesType k2CurvaturesDir(nPnts);
  ValuesType k1CurvaturesAng(nPnts);
  ValuesType k2CurvaturesAng(nPnts);
  double tmp1,tmp2,k1,k2;
  GradType grad;

  std::cout << "Starting curvature computation" << std::endl;

  std::cout << "   Isotropic odf" << std::endl;
  for( unsigned i = 0;i<nPnts;++i)
  {
    grad = grads->ElementAt(i);
    double theta = acos(grad[2]);
    double phi   = atan2(grad[1],grad[0]); // atan2(y,x) = atan(y/x);
    
    odf_i.ComputeCurvatures(theta,phi,tmp1,tmp2,k1,k2);
    k1CurvaturesAng[i] = k1 * odf_i[0] * vcl_sqrt( 1.0 / (4*vnl_math::pi) );
    k2CurvaturesAng[i] = k2 * odf_i[0] * vcl_sqrt( 1.0 / (4*vnl_math::pi) );
    if (isnan(k1))
    {
      std::cout << "NANs encountered" << std::endl;
      std::cout << grad << " - " << k1 << " " << k2 << std::endl;
      return EXIT_FAILURE;
    }
    if (! (areEqual(k1,k2) && areEqual(k1, (4*vnl_math::pi) ) ) )
    {
      std::cout << "Error curvatures are incorrect" << std::endl;
      std::cout << grad << " - " << k1 << " " << k2 << std::endl;
      return EXIT_FAILURE;      
    }
    std::cout << grad << " - " << k1 << " " << k2 << std::endl;

  }

  std::cout << "*************************************************************************************" << std::endl
            << "***  testCurvatureComputation NOT AN ACTUAL TEST NO COMPARISION MADE TO RESULTS  *****************************" << std::endl
            << "*************************************************************************************" << std::endl;

  return EXIT_FAILURE;
}

}// end EmptyNamespace
using namespace itk;
int itkSymRealSphericalHarmonicRepTest( int , char ** )
//int main( int argc, char * argv[] )
{
  if ( testNormalize() )
    return EXIT_FAILURE;
  if ( testEvaluation() )
    return EXIT_FAILURE;
  if ( testIndexing() )
    return EXIT_FAILURE;
  if ( testBasisFunctions() )
    return EXIT_FAILURE;
  if ( testCurvatureComputation() )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
