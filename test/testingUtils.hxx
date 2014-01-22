

#ifndef __testingUtils_hxx
#define __testingUtils_hxx

namespace DiffusionImagingTK_testing
{

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
