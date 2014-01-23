

#ifndef __testingUtils_hxx
#define __testingUtils_hxx

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMesh.h>
#include <itkRegularSphereMeshSource.h>
#include <iostream>
#include <fstream>
#include <itkImageRegionIteratorWithIndex.h>

namespace DiffusionImagingTK_testing
{

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
    PointType  point;

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

template <typename GradientDirectionContainerType>
void loadGradDirs( typename GradientDirectionContainerType::Pointer gradDirs, std::string bvecFile, std::string bvalFile )
{
  //Parse bvalFile
  std::vector<double> bValues;

  std::ifstream bvalIn(bvalFile.c_str());
  //first count the number of grad directions and make sure the only bvalues in the bval file
  // are 0 and a number...
  std::string line;
  while (! bvalIn.eof() )
  {
    getline(bvalIn,line);
    double val;
    std::stringstream ss(line);
    while (ss >> val)
    {
      bValues.push_back(val);
    }
  }
  
  //ok so now lets process the gradient table...
  //read in each line and put it in a string stream to process...
  std::ifstream bvecIn(bvecFile.c_str());
  getline(bvecIn,line);
  std::stringstream Xss(line);
  getline(bvecIn,line);
  std::stringstream Yss(line);
  getline(bvecIn,line);
  std::stringstream Zss(line);

  typename GradientDirectionContainerType::Element vect3d;
  unsigned int counter = 0;
  double x,y,z;
  double scale;
  while (Xss >> x)
  {
    scale = vcl_sqrt(bValues[counter]);
    Yss >> y;
    Zss >> z;
    vect3d[0] = scale * x; vect3d[1] = scale * y; vect3d[2] = scale * z;
    gradDirs->InsertElement( counter, vect3d );
    ++counter;
  }
  if (counter != bValues.size())
  {
    std::cerr << "different number of bvalues and gradients" << std::endl;
    gradDirs->Initialize();
  }
}

template < class DwiImageType4D, class DwiImageType >
void convertDWIimage(typename DwiImageType4D::Pointer inp4d, typename DwiImageType::Pointer outputIm )
{

  typedef itk::ImageRegionIteratorWithIndex< DwiImageType >         IteratorType;

  //Set up the gradient image size
  typename DwiImageType::SizeType  sizeGradImage;
  typename DwiImageType4D::SizeType size4D = inp4d->GetLargestPossibleRegion().GetSize();
  sizeGradImage[0] = size4D[0];
  sizeGradImage[1] = size4D[1];
  sizeGradImage[2] = size4D[2];
  outputIm->SetVectorLength(size4D[3]);

  typename DwiImageType::IndexType   indexGradImage = {{ 0, 0, 0 }};
  typename DwiImageType::RegionType  regionGradImage;
  regionGradImage.SetSize(  sizeGradImage );
  regionGradImage.SetIndex( indexGradImage);
  outputIm->SetRegions( regionGradImage );
  typename DwiImageType4D::SpacingType img4Dspacing = inp4d->GetSpacing();
  typename DwiImageType4D::PointType img4Dorigin = inp4d->GetOrigin();
  typename DwiImageType4D::DirectionType img4Ddir = inp4d->GetDirection();

  typename DwiImageType::SpacingType gradSpacing;
  typename DwiImageType::PointType gradOrigin;
  typename DwiImageType::DirectionType gradDirs;

  gradSpacing[0]  = img4Dspacing[0];  gradSpacing[1] = img4Dspacing[1];   gradSpacing[2] = img4Dspacing[2];
  gradOrigin[0]   = img4Dorigin[0];   gradOrigin[1]  = img4Dorigin[1];    gradOrigin[2] = img4Dorigin[2];

  for (unsigned int i = 0; i<3; ++i)
  {
    for (unsigned int j = 0; j<3; ++j)
    {
      gradDirs[i][j] = img4Ddir[i][j];
    }
  }

  outputIm->SetSpacing( gradSpacing );
  outputIm->SetOrigin( gradOrigin );
  outputIm->SetDirection( gradDirs );

  outputIm->Allocate();

  printf("Done GradIm->Allocate\n");

  //Copy data from img4d to gradim
  //THIS IS SLOW!!!
  IteratorType it( outputIm, outputIm->GetRequestedRegion() );

  typename DwiImageType::IndexType    gradIndex;
  typename DwiImageType::PixelType    gradPix;
  typename DwiImageType4D::IndexType  img4dIndex;

  //Probably a better way to do this but I don't really know what it is.
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    gradIndex = it.GetIndex();
    gradPix = it.Get();
    img4dIndex[0] = gradIndex[0];
    img4dIndex[1] = gradIndex[1];
    img4dIndex[2] = gradIndex[2];

    for ( unsigned int i=0; i<size4D[3]; ++i )
    {
      img4dIndex[3] = i;
      gradPix.SetElement( i, inp4d->GetPixel( img4dIndex ) );
    }
    it.Set( gradPix );
  }
  printf("Done Conversion\n");
}


}//End namespace

#endif
