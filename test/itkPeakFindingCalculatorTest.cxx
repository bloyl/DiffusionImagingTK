/**
 * @file  itkPeakFindingCalculatorTest.cxx
 * @brief Test itkPeakFindingCalculator module.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <itkSymRealSphericalHarmonicRep.h>
#include <itkPeakFindingCalculator.h>
#include <stdio.h>

#include <itkTestingMacros.h>

//#define  INFO_TO_FILE    // write all peaks to an ASCII file
#define  TEST_NALPLUT    // accelerated peak finding

const int    shrMaxOrder = 8;
const int    numberVoxls = 3;
const int    nMRtrxPeaks = 5;
const int    mrtrxBufSiz = nMRtrxPeaks * 3;
const int    mrtrxExtSiz = nMRtrxPeaks * 4;
const int    numShCoeffs = ( shrMaxOrder + 1 ) * ( shrMaxOrder + 2 ) / 2;
const double err_epsilon = 0.000002;

typedef itk::SymRealSphericalHarmonicRep < double, shrMaxOrder >  PixelType;


int testPeakFinding()
{
  // Luke SH --- Luke data
  double  shRepCoeffs[numberVoxls][numShCoeffs] =
          {
            {
               2.82094806e-01,  4.59387809e-01,  7.62639241e-03, -2.65089244e-01,
              -1.70459377e-03, -2.27291649e-03,  3.38489264e-01,  8.24082829e-03,
              -2.55495191e-01, -8.35784525e-03,  1.70923799e-01,  1.99015043e-03,
               1.02482527e-03, -5.59562538e-03, -3.80363222e-03,  1.70554623e-01,
               5.64062875e-03, -1.26915440e-01, -5.26009873e-03,  1.14953570e-01,
               4.45495034e-03, -7.90442452e-02, -1.32022612e-03, -3.45894368e-04,
               3.95073835e-03,  9.47112043e-04, -6.37689978e-03, -3.68007110e-03,
               4.89506535e-02,  2.26373016e-03, -3.73514183e-02, -2.18456332e-03,
               3.32550853e-02,  1.17909722e-03, -3.10148522e-02, -6.25214074e-04,
               2.12846715e-02,  3.60342383e-04,  3.37901074e-05, -1.33388466e-03,
              -1.39679192e-04,  2.47698487e-03,  2.21211842e-04, -3.72869126e-03,
              -2.16668053e-03
            },

            {
               2.82094806e-01, -4.56637770e-01, -2.07055407e-03, -2.66447693e-01,
              -3.83342733e-03, -2.03245622e-03,  3.29349637e-01,  6.67131925e-03,
               2.56829888e-01,  2.59662815e-03,  1.74307019e-01,  4.18567704e-03,
               2.51308928e-04,  3.66590638e-03,  3.40090320e-03, -1.61080286e-01,
              -6.84451079e-03, -1.24465309e-01, -5.27192699e-03, -1.18184127e-01,
              -1.89521688e-03, -8.28848481e-02, -2.20265985e-03,  4.79890703e-04,
              -2.14970903e-03,  4.34971036e-04, -1.83163770e-03, -3.47044342e-03,
               4.55690846e-02,  3.02502280e-03,  3.43799517e-02,  3.41011304e-03,
               3.30946259e-02,  2.10025162e-03,  3.35046090e-02,  6.53992000e-04,
               2.37989724e-02,  3.83429928e-04, -1.72163273e-04,  3.30030743e-04,
              -7.50337844e-04,  3.65657092e-04, -9.34689073e-04,  3.83923936e-04,
               2.20232457e-03
            },

            {
               2.82094806e-01, -1.84499600e-03, -4.37358674e-03, -2.57216394e-01,
              -2.21597566e-03, -5.69619704e-03,  3.87160391e-01, -5.81498304e-03,
              -4.88071050e-03,  3.55852279e-03,  1.55444667e-01,  1.25033641e-03,
               2.96612363e-03,  2.67166062e-03, -6.83821831e-03, -1.60360534e-03,
              -3.04495008e-03, -1.28530160e-01,  2.44779466e-03,  5.92144812e-03,
              -4.25606035e-04, -6.48434013e-02,  5.08837402e-04, -1.10527303e-03,
              -2.69610246e-05,  2.88396259e-03, -2.25038873e-03, -7.32504530e-03,
               1.14771903e-01, -3.58927739e-03, -2.93563283e-03,  1.99614253e-04,
               2.79063080e-02,  2.29189231e-04, -2.89570610e-03, -1.19774416e-03,
               1.54047897e-02, -7.50655425e-04,  2.04480690e-04, -1.01385545e-03,
              -9.12652526e-04, -5.92372544e-06,  1.25703122e-03,  2.75207148e-03,
              -3.13866278e-03
            }
          };


  #ifdef  TEST_NALPLUT         // Accelerated NALP computation
  double  allRawPeaks[numberVoxls][36] =
          {
            {
               0.999975,    -0.003092,     0.006333,    1.152481,
               0.492666,    -0.870218,     0.000241,    0.003035,
               0.503518,    -0.329527,     0.798675,    0.003025,
              -0.005536,    -0.562744,     0.826613,    0.002556,
               0.497721,     0.311159,     0.809601,    0.002063,
              -0.497486,    -0.867463,     0.004030,    0.002057,
              -0.007518,     0.543908,     0.839111,    0.001852,
              -0.016368,    -0.999862,     0.002858,    0.001511,
               0.004595,     0.001553,     0.999988,    0.001381
            },

            {
               0.002523,    -0.999994,    -0.002558,    1.139916,
               0.202992,    -0.526240,     0.825752,    0.004237,
               0.266560,     0.513826,     0.815432,    0.003746,
              -0.859480,     0.511146,    -0.004907,    0.002782,
               0.003691,     0.000778,     0.999993,    0.002400,
               0.552324,    -0.005179,     0.833613,    0.002388,
               0.999806,     0.019642,     0.001253,    0.001850,
              -0.535278,     0.013071,     0.844575,    0.001308,
               0.854526,     0.519402,     0.002704,    0.001117
            },

            {
               0.999953,    -0.008272,    -0.005101,    0.645980,
               0.000102,    -0.999998,    -0.001953,    0.633885,
               0.496562,     0.014951,     0.867873,    0.001673,
              -0.019403,    -0.473921,     0.880354,    0.001405,
              -0.026013,     0.012218,     0.999587,    0.000969,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000
            }
          };
  #else                        // brute-force NALP computation
  double  allRawPeaks[numberVoxls][36] =
          {
            {
               0.999975,    -0.003092,     0.006333,    1.152545,
               0.492664,    -0.870220,     0.000242,    0.003035,
               0.503522,    -0.329606,     0.798640,    0.003026,
              -0.005536,    -0.562743,     0.826614,    0.002556,
               0.497724,     0.311133,     0.809609,    0.002064,
              -0.497484,    -0.867464,     0.004030,    0.002057,
              -0.007519,     0.543919,     0.839104,    0.001852,
              -0.016368,    -0.999862,     0.002859,    0.001511,
               0.004601,     0.000860,     0.999989,    0.001381
            },

            {
               0.002523,    -0.999994,    -0.002558,    1.139935,
               0.202934,    -0.526256,     0.825756,    0.004239,
               0.266477,     0.513847,     0.815445,    0.003748,
              -0.859481,     0.511145,    -0.004907,    0.002782,
               0.002925,     0.001083,     0.999995,    0.002400,
               0.552332,    -0.005179,     0.833608,    0.002388,
               0.999806,     0.019643,     0.001253,    0.001850,
              -0.535287,     0.013072,     0.844569,    0.001308,
               0.854526,     0.519401,     0.002704,    0.001117
            },

            {
               0.999953,    -0.008272,    -0.005100,    0.646009,
               0.000102,    -0.999998,    -0.001953,    0.633904,
               0.496574,     0.014955,     0.867865,    0.001674,
              -0.019408,    -0.473928,     0.880350,    0.001406,
              -0.026019,     0.012180,     0.999587,    0.000969,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000,
               0.000000,     0.000000,     0.000000,    0.000000
            }
          };
  #endif


  #ifdef INFO_TO_FILE
  // MRtrix peaks --- Luke data (simulation)
  // voxel #0: horizontal --- [ 1 0 0 ]
  // voxel #1: vertical   --- [ 0 1 0 ]
  // voxel #2: cross      --- [ 1 0 0 ] + [ 0 1 0 ]
  double  mrtrixPeaks[numberVoxls][mrtrxBufSiz] =
          {
            {
               7.00846732e-01, -2.17178208e-03,  4.20162128e-03,
               6.21720329e-02, -1.09789588e-01,  7.60147523e-05,
              -6.25087544e-02, -1.09033048e-01,  3.15546022e-05,
              -2.04979209e-03, -1.25389695e-01,  7.78248213e-05,
              -1.04265660e-03, -3.31950927e-04,  1.71643763e-03
            },

            {
               1.75462803e-03, -6.96365118e-01, -1.70552416e-03,
              -1.09852873e-01,  6.53061271e-02, -1.26035549e-04,
               1.27315268e-01,  2.50503793e-03, -1.87696187e-05,
               1.08514570e-01,  6.59241900e-02,  1.20940895e-04,
               5.00256312e-04, -1.28934369e-03,  2.02249736e-03
            },

            {
               4.40724611e-01, -3.64602846e-03, -1.94008509e-03,
               4.86247700e-05, -4.34704721e-01, -7.07804633e-04,
              -8.33624750e-02, -8.26655477e-02, -2.46377051e-04,
              -2.39618821e-05,  1.13826654e-05,  8.88229639e-04,
               3.83511477e-04,  1.08955828e-05,  6.58161531e-04
            }
          };
  #endif


  #ifdef INFO_TO_FILE
  FILE   *   outFile = fopen( "AllPeaks.txt", "w" );
  #endif


  int        i, k;
  std::vector< double >            odPeaks; // each peak by 4 components: vx, vy, vz, mag
  PixelType::GradientDirectionType peakDir;

  typedef itk::PeakFindingCalculator < PixelType, double, double > PeakFindingCalculator;
  PeakFindingCalculator::Pointer       pFinder = PeakFindingCalculator::New();
  EXERCISE_BASIC_OBJECT_METHODS( pFinder, PeakFindingCalculator );

  PeakFindingCalculator::InputType     shCoefs[ numberVoxls ];
  for ( i = 0; i < numberVoxls; i ++ ) shCoefs[i] = shRepCoeffs[i];


  for ( k = 0; k < numberVoxls; k ++ )
  {
    odPeaks.clear();

    // ------------------------------- the test -------------------------------

    pFinder->SetCoefficients( shCoefs[k] );
    pFinder->SetCoefficient( 0, shCoefs[k][0] );

    pFinder->SetPeakFiltering( 0 );
    #ifdef TEST_NALPLUT
    pFinder->SetAcceleration( 1 );
    #else
    pFinder->SetAcceleration( 0 );
    #endif
    int   nmPeaks   = pFinder->GetRawPeaks( odPeaks );
    nmPeaks += 0; // just to suppress warnings (var only used for INFO_TO_FILE)


    int   numVals   = (  ( k == 2 ) ? ( 5 * 4 ) : ( 9 * 4 )  );
    if (  numVals  !=  int( odPeaks.size() )  )
    {
      printf( "\nError with the number of peaks for voxel #%d\n", k );
      #ifdef INFO_TO_FILE
      fclose( outFile );
      outFile = NULL;
      #endif
      return 1;
    }

    // zero padding for voxel #2
    for ( i = numVals; i < 36; i ++ ) odPeaks.push_back( 0.0 );

    for ( i = 0; i < 36; i ++ )
    {
      if (  vcl_fabs( odPeaks[i] - allRawPeaks[k][i] )  >  err_epsilon  )
      {
        printf( "\nError with the peak result of voxel #%d, i = %d\n", k, i );
        #ifdef  INFO_TO_FILE
        fclose( outFile );
        outFile = NULL;
        #endif
        return 1;
      }
    }


    // -------------------------------------------------------------------------


    #ifdef   INFO_TO_FILE
    fprintf( outFile, "\n\n================ Voxel %d ================\n\n", k );

    // show the ground truth (for Luke simulation data)
    peakDir.fill( 0 );
    PixelType odShRep( shRepCoeffs[k] );
    fprintf ( outFile, "-------- ground truth --------\n" );

    if ( k == 0 || k == 2 )
    {
      peakDir[0] = 1;
      fprintf(  outFile,  "%+.6f    %+.6f    %+.6f  =  %+.6f\n",
                peakDir[0],  peakDir[1],  peakDir[2],
                odShRep.Evaluate( peakDir )  );
      peakDir.fill( 0 );
    }

    if ( k == 1 || k == 2 )
    {
      peakDir[1] = 1;
      fprintf(  outFile,  "%+.6f    %+.6f    %+.6f  =  %+.6f\n",
                peakDir[0],  peakDir[1],  peakDir[2],
                odShRep.Evaluate( peakDir )  );
      peakDir.fill( 0 );
    }

    // show the extracted peaks
    fprintf( outFile, "\n-------- estimated peaks --------\n");
    for ( i = 0; i < nmPeaks; i ++ )
    {
      int  index = i * 4;
      peakDir[0] = odPeaks[ index     ];
      peakDir[1] = odPeaks[ index + 1 ];
      peakDir[2] = odPeaks[ index + 2 ];
      fprintf(  outFile,  "%+.6f    %+.6f    %+.6f    %+.6f -------- %+.6f\n",
                peakDir[0],  peakDir[1],  peakDir[2],  odPeaks[ index + 3 ],
                odShRep.Evaluate( peakDir )  );
    }

    // use a new array to save 4-component tuples: vx, vy, vz, vecMag
    double  MRtrixPeaks[ mrtrxExtSiz ];
    for ( j = 0; j < nMRtrxPeaks; j ++ )
    {
      int    index0 = j * 3;
      int    index1 = j * 4;
      double vecMag =
             vcl_sqrt(  mrtrixPeaks[k][ index0     ] * mrtrixPeaks[k][ index0     ] +
                        mrtrixPeaks[k][ index0 + 1 ] * mrtrixPeaks[k][ index0 + 1 ] +
                        mrtrixPeaks[k][ index0 + 2 ] * mrtrixPeaks[k][ index0 + 2 ]  );
      double gamCev = ( vecMag == 0.0 ) ? 1.0: ( 1.0 / vecMag );

      for ( i = 0; i < 3; i ++ )
      MRtrixPeaks[ index1 + i ] = mrtrixPeaks[k][ index0 + i ] * gamCev;

      MRtrixPeaks[ index1 + 3 ] = vecMag;
    }

    // sort the MRtrix result by the OD value in decreasing order
    for ( j = 0; j < nMRtrxPeaks - 1; j ++ )
    {
      int index0 = j * 4;
      for ( i = j + 1; i < nMRtrxPeaks; i ++ )
      {
        int index1 = i * 4;
        if (  MRtrixPeaks[ index0 + 3 ]  <  MRtrixPeaks[ index1 + 3 ]  )
        {
          for ( int n = 0; n < 4; n ++ )
          {
            double  tmpVal = MRtrixPeaks[ index0 + n ];
            MRtrixPeaks[ index0 + n ] = MRtrixPeaks[ index1 + n ];
            MRtrixPeaks[ index1 + n ] = tmpVal;
          }
        }
      }
    }

    // show the MRtrix result
    fprintf( outFile, "\n-------- mrtrix peaks --------\n" );
    for ( i = 0; i < mrtrxExtSiz; i += 4 )
    {
      peakDir[0] = MRtrixPeaks[ i     ];
      peakDir[1] = MRtrixPeaks[ i + 1 ];
      peakDir[2] = MRtrixPeaks[ i + 2 ];
      fprintf(  outFile,  "%+.6f    %+.6f    %+.6f    %+.6f -------- %+.6f\n",
                peakDir[0],  peakDir[1],  peakDir[2],  MRtrixPeaks[ i + 3 ],
                odShRep.Evaluate( peakDir )  );
    }
    #endif

  } // endfor: each voxel

  #ifdef  INFO_TO_FILE
  fclose( outFile );
  outFile = NULL;
  #endif

  return 0;
}

using namespace itk;
int itkPeakFindingCalculatorTest( int , char ** )
{
  if ( testPeakFinding() ) return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
