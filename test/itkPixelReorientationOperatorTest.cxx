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

#include <iostream>
#include <iomanip>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_matrix_fixed.h>
#include <itkArray2D.h>
#include <itkMatrix.h>
#include <vnl/vnl_random.h>
#include <itkTimeProbe.h>
#include <itkTestingMacros.h>

#include <itkPixelReorientationOperator.h>
#include <itkDTIPixelReorientationOperator.h>
#include <itkRSHPixelReorientationOperator.h>

const double percision = 0.0001;
const unsigned int numTests = 100;

vnl_matrix<double>
mkRotMat(double a, double b, double g)
{
    vnl_matrix<double> rotA(3,3,0);
    vnl_matrix<double> rotB(3,3,0);
    vnl_matrix<double> rotG(3,3,0);

    rotA(0,0) = vcl_cos(a); rotA(0,1) = -vcl_sin(a);
    rotA(1,0) = vcl_sin(a); rotA(1,1) = vcl_cos(a);
    rotA(2,2) = 1;

    rotG(0,0) = vcl_cos(g); rotG(0,1) = -vcl_sin(g);
    rotG(1,0) = vcl_sin(g); rotG(1,1) = vcl_cos(g);
    rotG(2,2) = 1;

    rotB(0,0) = vcl_cos(b); rotB(0,2) = vcl_sin(b);
    rotB(1,1) = 1;
    rotB(2,0) = -vcl_sin(b); rotB(2,2) = vcl_cos(b);

    vnl_matrix<double> res = rotG * rotB * rotA;
    return res;
}

template <typename TDtiType, typename TDirType>
double evaluateDTI(TDtiType dti, TDirType dir)
{
  double val = 0;
  val +=     dti[0] * dir[0] * dir[0];
  val += 2 * dti[1] * dir[0] * dir[1];
  val += 2 * dti[2] * dir[0] * dir[2];
  val +=     dti[3] * dir[1] * dir[1];
  val += 2 * dti[4] * dir[1] * dir[2];
  val +=     dti[5] * dir[2] * dir[2];
  return val;
}

int itkPixelReorientationOperationOnScalarsTest( vnl_random )
{
  typedef float                                        PixelType;
  typedef  vnl_matrix<double>                          JacobianType;

  typedef itk::PixelReorientationOperator<PixelType>   OrienterType;

  OrienterType::Pointer orienter  = OrienterType::New();
  EXERCISE_BASIC_OBJECT_METHODS( orienter, OrienterType );
  
  JacobianType rotMat = mkRotMat(0.1,0.12,0.23);

  PixelType pix = 1.0;
  PixelType pix2 = orienter->Reorient( pix, rotMat );

  if (pix != pix2)
  {
    std::cout << "Failed Pixel Reorientation test on scalars" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "itkPixelReorientationOperationOnScalarsTest Passed\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnDTITest( vnl_random )
{
  typedef itk::DiffusionTensor3D<float>                 PixelType;

  typedef itk::PixelReorientationOperator<PixelType>    OrienterType;

  OrienterType::Pointer orienter  = OrienterType::New();
  EXERCISE_BASIC_OBJECT_METHODS( orienter, OrienterType );

  vnl_matrix<double> rotMat = mkRotMat(0.1,0.12,0.23);

  PixelType inputPixel;
  PixelType outputPixel;

  inputPixel(0,0) =  19.0;
  inputPixel(0,1) =   0.0;
  inputPixel(0,2) =   0.0;
  inputPixel(1,0) =   0.0; // overrides (0,1)
  inputPixel(1,1) =  23.0;
  inputPixel(1,2) =   0.0;
  inputPixel(2,0) =   7.0; // overrides (0,2)
  inputPixel(2,1) =   0.0; // overrides (1,2)
  inputPixel(2,2) =  29.0;

  // Correct rotation of tensor3D by matrix3D
  outputPixel[0] =  21.086500;
  outputPixel[1] =  -0.749525;
  outputPixel[2] =   7.570830;
  outputPixel[3] =  22.711900;
  outputPixel[4] =   2.534080;
  outputPixel[5] =  27.201700;

  PixelType pix2 = orienter->Reorient( inputPixel, rotMat );

  for (unsigned int i=0; i<6; ++i)
  {
    if ( vcl_abs(pix2[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on Dti" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //test reorientation with Array2d
  itk::Array2D<double> arrayMat = static_cast< itk::Array2D<double> >(rotMat);
  PixelType pix3 = orienter->Reorient( inputPixel, arrayMat );
  for (unsigned int i=0; i<6; ++i)
  {
    if ( vcl_abs(pix3[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on Dti" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //test reorientation with itkMatrix
  itk::Matrix<double> itkMat = static_cast< itk::Matrix<double> >(rotMat);
  PixelType pix4 = orienter->Reorient( inputPixel, itkMat );
  for (unsigned int i=0; i<6; ++i)
  {
    if ( vcl_abs(pix4[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on Dti" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //Add failure test for non 3x3 mats

  std::cout << "itkPixelReorientationOperationOnDTITest Passed\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnDTITest2( vnl_random rand )
{

  typedef itk::DiffusionTensor3D<float>                 PixelType;

  typedef itk::PixelReorientationOperator<PixelType>    OrienterType;

  OrienterType::Pointer orienter  = OrienterType::New();

  PixelType pix1;
  pix1[0] = rand.drand32(-100,100);
  pix1[1] = rand.drand32(-100,100);
  pix1[2] = rand.drand32(-100,100);
  pix1[3] = rand.drand32(-100,100);
  pix1[4] = rand.drand32(-100,100);
  pix1[5] = rand.drand32(-100,100);

  typedef vnl_vector_fixed< double, 3 >                GradType;

  vnl_matrix<double> rotMat;

  PixelType rotPix1;

  GradType dir;

  double alpha,beta,gamma;
  double sign;

  itk::TimeProbe timer;
  timer.Start();

  double maxError = 0;
  double error = 0;
  for (unsigned int i = 0;i < numTests;++i)
  {

    dir[0] = rand.drand32(-1,1);
    dir[1] = rand.drand32(-1,1);
    dir[2] = rand.drand32(-1,1);
    dir = dir/dir.two_norm();

    //Pick these at random!
    alpha = rand.drand32(-3.14159,3.14159);
    beta = rand.drand32(-1.5708,1.5708);
    gamma = rand.drand32(-3.14159,3.14159);
    sign = rand.drand32(-1,1);

    rotMat = mkRotMat(alpha,beta,gamma);
    //add reflections
    if (sign < 0)
    {
      rotMat *= -1;
    }

    double val1 = evaluateDTI(pix1,rotMat.transpose() * dir);
    error = vcl_abs( (val1 - evaluateDTI(orienter->Reorient(pix1,rotMat),dir)) /val1);

    if (error > maxError) maxError = error;
    if ( error > 0.0001 )
    {
      std::cout << "Rotation Test Failed!" << std::endl;
      std::cout << "i             :" << i << std::endl;
      std::cout << "rotMat        :\n" << rotMat << std::endl;
      std::cout << "dir           : " << dir << std::endl;
      std::cout << "rotdir        : " << rotMat * dir << std::endl;

      std::cout << "\nalpha           : " << alpha << std::endl;
      std::cout << "gamma           : " << gamma << std::endl;

      std::cout << "evalutate     : " << evaluateDTI(pix1,rotMat.transpose() * dir) << std::endl;
      std::cout << "Rot Evaluate  : " << evaluateDTI(orienter->Reorient(pix1,rotMat),dir) << std::endl;
      std::cout << "error         : " << error << std::endl;

      std::cout << "itkPixelReorientationOperationOnDTITest2 (Reorient)  FAILED\n";
      return EXIT_FAILURE;
    }
  }
  timer.Stop();
  std::cout << "RotationTest: Time Elapsed " << timer.GetMean() << " seconds.\n";
  std::cout << "RotationTest: Time per Rotation " << timer.GetMean() / numTests << " seconds.\n";
  std::cout << "max percent Error : " << maxError << "\n";
  std::cout << "itkPixelReorientationOperationOnDTITest2 (Reorient)  PASSED\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnDTITest3( vnl_random rand )
{

  typedef itk::DiffusionTensor3D<float>                 PixelType;

  typedef itk::PixelReorientationOperator<PixelType>    OrienterType;

  OrienterType::Pointer orienter  = OrienterType::New();

  PixelType pix1;

  pix1[0] = rand.drand32(-100,100);
  pix1[1] = rand.drand32(-100,100);
  pix1[2] = rand.drand32(-100,100);
  pix1[3] = rand.drand32(-100,100);
  pix1[4] = rand.drand32(-100,100);
  pix1[5] = rand.drand32(-100,100);

  typedef vnl_vector_fixed< double, 3 >                GradType;

  vnl_matrix<double> rotMat;

  PixelType rotPix1;

  GradType dir;

  double alpha,beta,gamma;

  itk::TimeProbe timer;
  timer.Start();

  double maxError = 0;
  double error = 0;
  for (unsigned int i = 0;i < numTests;++i)
  {

    dir[0] = rand.drand32(-1,1);
    dir[1] = rand.drand32(-1,1);
    dir[2] = rand.drand32(-1,1);
    dir = dir/dir.two_norm();

    //Pick these at random!
    alpha = rand.drand32(-3.14159,3.14159);
    beta = rand.drand32(-1.5708,1.5708);
    gamma = rand.drand32(-3.14159,3.14159);

    rotMat = mkRotMat(alpha,beta,gamma);
    double val1 = evaluateDTI(pix1,rotMat * dir);
    error = vcl_abs( (val1 - evaluateDTI(orienter->ReorientInv(pix1,rotMat),dir)) /val1);

    if (error > maxError) maxError = error;
    if ( error > 0.0001 )
    {
      std::cout << "Rotation Test Failed!" << std::endl;
      std::cout << "i             :" << i << std::endl;
      std::cout << "rotMat        :\n" << rotMat << std::endl;
      std::cout << "dir           : " << dir << std::endl;
      std::cout << "rotdir        : " << rotMat * dir << std::endl;

      std::cout << "\nalpha           : " << alpha << std::endl;
      std::cout << "gamma           : " << gamma << std::endl;

      std::cout << "evalutate     : " << evaluateDTI(pix1,rotMat * dir) << std::endl;
      std::cout << "Rot Evaluate  : " << evaluateDTI(orienter->ReorientInv(pix1,rotMat),dir) << std::endl;
      std::cout << "error         : " << error << std::endl;

      std::cout << "itkPixelReorientationOperationOnDTITest3 (ReorientInv)  FAILED\n";
      return EXIT_FAILURE;
    }
  }
  timer.Stop();
  std::cout << "RotationTest: Time Elapsed " << timer.GetMean() << " seconds.\n";
  std::cout << "RotationTest: Time per Rotation " << timer.GetMean() / numTests << " seconds.\n";
  std::cout << "max percent Error : " << maxError << "\n";
  std::cout << "itkPixelReorientationOperationOnDTITest3 (ReorientInv)  PASSED\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnRSHTest( vnl_random )
{
  typedef itk::SymRealSphericalHarmonicRep<float,8>   PixelType;
  typedef vnl_matrix<double>                          JacobianType;

  typedef itk::PixelReorientationOperator<PixelType>  OrienterType;
  OrienterType::Pointer orienter  = OrienterType::New();
  EXERCISE_BASIC_OBJECT_METHODS( orienter, OrienterType );

  JacobianType rotMat = mkRotMat(0.1,0.12,0.23);

  float c1[45] = {0.7212520319878115,-0.0003385282051116059,0.0000749605643584339,-0.595677776291134,0.00030204505569102896,
   0.0007960419344463232,0.0007041597783487286,-0.00004763110564954577,0.00023051855355286227,-0.00041468059636469434,
   0.31454381107471285,0.0003240763081598791,-0.0003646591186538331,-0.00013836214560059706,0.0010292282395887864,
   0.0011945036277760113,0.00044623369448974076,-0.0007192376368217943,0.0001093875566357573,-0.00015953304798718704,
   -0.0014569479866522461,-0.12532455607008042,0.001029082609176401,0.00015452420774457638,-0.0002648896339109974,
   -0.0012802266647137153,0.000852198927608263,0.001338962813971084,0.003603085155208449,-0.00008088348021568353,
   -0.0006703905742253685,0.000771066948178821,0.0012535994711646709,-0.0012808114574620819,0.000998683995279085,
   0.0014864901181477208,0.03818977726185059,-0.0015034416680695956,-0.0018179313371678475,-0.0020668893003337463,
   0.0019221161938821688,-0.0010735703142327364,-0.0013657370782731923,-0.00015421727639094215,-0.0003780481257301127
  };

  float c2[45] = {0.721252,  -0.00734784,  -0.119153,  -0.5829,  0.0282855,  -0.00288973,  -0.000797906,  0.000918557,
    0.0136837,  0.111875 , 0.292442 , -0.0259144 , 0.00634702 , -0.00053799 , 0.000962706 , -0.00149031 , 0.00140451,
    0.000528821,  -0.00160843, -0.011525 , -0.0629438 , -0.106401 , 0.0158689 , -0.00626362 , 0.000496838 , -0.00102716,
    -8.70135e-05,  0.000677845,  -0.00285198,  0.000876164,  0.000646569,  -0.00239911,  -0.00148314,  0.000165997,
    0.00826977,  0.0230242,  0.0283441,  -0.00777702,  0.00306007,  -8.53282e-05,  0.000888803,  -0.00113396,  0.000301292,
    0.000684019 , 0.0019872};

  PixelType inputPixel(c1);
  PixelType outputPixel(c2);

  //~ std::cout << "Expected output : "<< outputPixel << std::endl;

  PixelType pix2 = orienter->Reorient( inputPixel, rotMat );
  //~ std::cout << "vnlMatrix rotMat output: "<< pix2 << std::endl;
  for (unsigned int i=0; i<45; ++i)
  {
    if ( vcl_abs(pix2[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on RSH" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //test reorientation with Array2d
  itk::Array2D<double> arrayMat = static_cast< itk::Array2D<double> >(rotMat);
  PixelType pix3 = orienter->Reorient( inputPixel, arrayMat );
  //~ std::cout << "array2d rotMat output: "<< pix3 << std::endl;
  for (unsigned int i=0; i<45; ++i)
  {
    if ( vcl_abs(pix3[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on scalars" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //test reorientation with Array2d
  itk::Matrix<double> itkMat = static_cast< itk::Matrix<double> >(rotMat);
  PixelType pix4 = orienter->Reorient( inputPixel, itkMat );
  //~ std::cout << "itk matrix rotMat output: "<< pix4 << std::endl;
  for (unsigned int i=0; i<45; ++i)
  {
    if ( vcl_abs(pix4[i] - outputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on scalars" << std::endl;
      return EXIT_FAILURE;
    }
  }

  //Add failure test for non 3x3 mats

  //Test special Cases of rotMat
  JacobianType rotId = mkRotMat(0.0,0.0,0.0);
  PixelType pix5 = orienter->Reorient( inputPixel, rotId );
  //~ std::cout << "vnlMatrix rotMat output: "<< pix2 << std::endl;
  for (unsigned int i=0; i<45; ++i)
  {
    if ( vcl_abs(pix5[i] - inputPixel[i]) > percision )
    {
      std::cout << "Failed Pixel Reorientation test on RSH - Identity Matrix" << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnRSHTest2( vnl_random rand )
{

  typedef double           PrecisionType;

  PrecisionType c1[45] = {0.7212520319878115,-0.0003385282051116059,0.0000749605643584339,-0.595677776291134,0.00030204505569102896,
   0.0007960419344463232,0.0007041597783487286,-0.00004763110564954577,0.00023051855355286227,-0.00041468059636469434,
   0.31454381107471285,0.0003240763081598791,-0.0003646591186538331,-0.00013836214560059706,0.0010292282395887864,
   0.0011945036277760113,0.00044623369448974076,-0.0007192376368217943,0.0001093875566357573,-0.00015953304798718704,
   -0.0014569479866522461,-0.12532455607008042,0.001029082609176401,0.00015452420774457638,-0.0002648896339109974,
   -0.0012802266647137153,0.000852198927608263,0.001338962813971084,0.003603085155208449,-0.00008088348021568353,
   -0.0006703905742253685,0.000771066948178821,0.0012535994711646709,-0.0012808114574620819,0.000998683995279085,
   0.0014864901181477208,0.03818977726185059,-0.0015034416680695956,-0.0018179313371678475,-0.0020668893003337463,
   0.0019221161938821688,-0.0010735703142327364,-0.0013657370782731923,-0.00015421727639094215,-0.0003780481257301127
  };

  typedef itk::SymRealSphericalHarmonicRep< PrecisionType, 8 >
                                                          PixelType;

  typedef itk::PixelReorientationOperator<PixelType>      OrienterType;
  OrienterType::Pointer orienter  = OrienterType::New();

  typedef PixelType::GradientDirectionType                GradType;

  vnl_matrix<double> rotMat = mkRotMat(0.1,0.12,0.23);

  PixelType pix1(c1);
  pix1 *= 1000;
  PixelType rotPix1;

  GradType dir;

  double alpha,beta,gamma;

  itk::TimeProbe timer;
  timer.Start();

  double maxError = 0;
  double error = 0;
  double tmp;
  for (unsigned int i = 0;i < numTests;++i)
  {
    dir[0] = rand.drand32(-1,1);
    dir[1] = rand.drand32(-1,1);
    dir[2] = rand.drand32(-1,1);
    dir = dir/dir.two_norm();

    //Pick these at random!
    alpha = rand.drand32(-3.14159,3.14159);
    gamma = rand.drand32(-3.14159,3.14159);

    //For 10% beta=0 and 10% beta=pi
    tmp = rand.drand32(0,1);
    if (tmp < 0.1)
    {
      beta = 0;
    }
    else if (tmp < 0.2)
    {
      beta = vnl_math::pi;
    }
    else
    {
      beta = rand.drand32(-1.5708,1.5708);
    }

    rotMat = mkRotMat(alpha,beta,gamma);
    double val1 = pix1.Evaluate(rotMat.transpose() * dir);
    error = vcl_abs( (val1 - (orienter->Reorient(pix1,rotMat)).Evaluate(dir)) /val1);
    if (error > maxError) maxError = error;
    if ( error > 0.0001 )
    {
      std::cout << "i             :" << i << std::endl;
      std::cout << "rotMat        :\n" << rotMat << std::endl;
      std::cout << "rotMat2       :\n" << mkRotMat(0.0,0.0,2.81422) << std::endl;
      std::cout << "dir           : " << dir << std::endl;
      std::cout << "rotdir        : " << rotMat * dir << std::endl << std::endl;


      std::cout << "alpha         : " << alpha << std::endl;
      std::cout << "beta          : " << beta << std::endl;
      std::cout << "gamma         : " << gamma << std::endl;

      // std::cout << "evalutate     : " << pix1.Evaluate(rotMat.transpose() * dir) << std::endl;
      // std::cout << "Rot Evaluate  : " << orienter->Reorient(pix1, rotMat).Evaluate(dir) << std::endl;
      std::cout << "error         : " << error << std::endl;
      std::cout << "itkPixelReorientationOperationOnRSHTest2  (Reorient)     FAILED\n";
      return EXIT_FAILURE;
    }
  }

  timer.Stop();
  std::cout << "RotationTest: Time Elapsed " << timer.GetMean() << " seconds.\n";
  std::cout << "RotationTest: Time per Rotation " << timer.GetMean() / numTests << " seconds.\n";
  std::cout << "max percent Error : " << maxError << "\n";
  std::cout << "itkPixelReorientationOperationOnRSHTest2  (Reorient)          PASSED\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperationOnRSHTest3( vnl_random rand )
{

  typedef double           PrecisionType;

  PrecisionType c1[45] = {0.7212520319878115,-0.0003385282051116059,0.0000749605643584339,-0.595677776291134,0.00030204505569102896,
   0.0007960419344463232,0.0007041597783487286,-0.00004763110564954577,0.00023051855355286227,-0.00041468059636469434,
   0.31454381107471285,0.0003240763081598791,-0.0003646591186538331,-0.00013836214560059706,0.0010292282395887864,
   0.0011945036277760113,0.00044623369448974076,-0.0007192376368217943,0.0001093875566357573,-0.00015953304798718704,
   -0.0014569479866522461,-0.12532455607008042,0.001029082609176401,0.00015452420774457638,-0.0002648896339109974,
   -0.0012802266647137153,0.000852198927608263,0.001338962813971084,0.003603085155208449,-0.00008088348021568353,
   -0.0006703905742253685,0.000771066948178821,0.0012535994711646709,-0.0012808114574620819,0.000998683995279085,
   0.0014864901181477208,0.03818977726185059,-0.0015034416680695956,-0.0018179313371678475,-0.0020668893003337463,
   0.0019221161938821688,-0.0010735703142327364,-0.0013657370782731923,-0.00015421727639094215,-0.0003780481257301127
  };

  typedef itk::SymRealSphericalHarmonicRep< PrecisionType, 8 >
                                                          PixelType;

  typedef itk::PixelReorientationOperator<PixelType>      OrienterType;
  OrienterType::Pointer orienter  = OrienterType::New();

  typedef PixelType::GradientDirectionType                GradType;

  vnl_matrix<double> rotMat = mkRotMat(0.1,0.12,0.23);

  PixelType pix1(c1);
  pix1 *= 1000;
  PixelType rotPix1;

  GradType dir;

  double alpha,beta,gamma;

  itk::TimeProbe timer;
  timer.Start();

  double maxError = 0;
  double error = 0;
  for (unsigned int i = 0;i < numTests;++i)
  {

    dir[0] = rand.drand32(-1,1);
    dir[1] = rand.drand32(-1,1);
    dir[2] = rand.drand32(-1,1);
    dir = dir/dir.two_norm();

    //Pick these at random!
    alpha = rand.drand32(-3.14159,3.14159);
    beta = rand.drand32(-1.5708,1.5708);
    gamma = rand.drand32(-3.14159,3.14159);

    rotMat = mkRotMat(alpha,beta,gamma);
    double val1 = pix1.Evaluate(rotMat * dir);
    error = vcl_abs( (val1 - (orienter->ReorientInv(pix1,rotMat)).Evaluate(dir)) /val1);
    if (error > maxError) maxError = error;
    if ( error > 0.0001 )
    {
      std::cout << "i             :" << i << std::endl;
      std::cout << "rotMat        :\n" << rotMat << std::endl;
      std::cout << "dir           : " << dir << std::endl;
      std::cout << "rotdir        : " << rotMat * dir << std::endl;

      std::cout << "\nalpha           : " << alpha << std::endl;
      std::cout << "gamma           : " << gamma << std::endl;

      std::cout << "evalutate     : " << pix1.Evaluate(rotMat * dir) << std::endl;
      std::cout << "Rot Evaluate  : " << orienter->ReorientInv(pix1, rotMat).Evaluate(dir) << std::endl;
      std::cout << "error         : " << error << std::endl;
      std::cout << "itkPixelReorientationOperationOnRSHTest3 (ReorientInv)      FAILED\n";
      return EXIT_FAILURE;
    }
  }

  timer.Stop();
  std::cout << "RotationTest: Time Elapsed " << timer.GetMean() << " seconds.\n";
  std::cout << "RotationTest: Time per Rotation " << timer.GetMean() / numTests << " seconds.\n";
  std::cout << "max percent Error : " << maxError << "\n";
  std::cout << "itkPixelReorientationOperationOnRSHTest3  (ReorientInv)       PASSED\n\n";
  return EXIT_SUCCESS;
}

int itkPixelReorientationOperatorTest( int , char ** )
{

  const long seed = 10000L;
  vnl_random rand(seed);

  if (itkPixelReorientationOperationOnScalarsTest(rand) )
    return EXIT_FAILURE;

  if (itkPixelReorientationOperationOnDTITest(rand) )
    return EXIT_FAILURE;

  if (itkPixelReorientationOperationOnDTITest2(rand) )
    return EXIT_FAILURE;

  if (itkPixelReorientationOperationOnDTITest3(rand) )
    return EXIT_FAILURE;

  if (itkPixelReorientationOperationOnRSHTest(rand) )
    return EXIT_FAILURE;

 if (itkPixelReorientationOperationOnRSHTest2(rand) )
   return EXIT_FAILURE;

  if (itkPixelReorientationOperationOnRSHTest3(rand) )
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}

