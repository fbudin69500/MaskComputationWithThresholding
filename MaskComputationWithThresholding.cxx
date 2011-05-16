#include "MaskComputationWithThresholdingCLP.h"
#include "itkImageFileWriter.h"
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkGrayscaleErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkCastImageFilter.h>
#include <itkScalarImageToHistogramGenerator.h>
#include <itkOtsuMultipleThresholdsCalculator.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <limits>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1



template<class T> int DoIt( int argc, char * argv[] )
{
    /////////////////   Get arguments using GenerateCLP parser to get all the arguments  ///////////////////////
    PARSE_ARGS ;
    typedef itk::Image< T, 3 > InputImageType ;
    /////////////////   Load input file    ////////////////////////////////////////////////
    typedef itk::ImageFileReader< InputImageType > ImageReaderType ;
    typename ImageReaderType::Pointer reader = ImageReaderType::New() ;
    reader->SetFileName( inputVolume ) ;
    reader->Update() ;
    typename InputImageType::Pointer image ;
    image = reader->GetOutput() ;
    ////////Automatic threshold////////////////////////////
    if( threshold || autoThreshold )
    {
      typedef itk::MinimumMaximumImageCalculator< InputImageType > MinMaxFilterType ;
      typename MinMaxFilterType::Pointer minMaxFilter = MinMaxFilterType::New() ;
      minMaxFilter->SetImage( image ) ;
      minMaxFilter->Compute() ;
      typedef itk::Statistics::ScalarImageToHistogramGenerator< InputImageType > ImageToHistogramType ;;
      typename ImageToHistogramType::Pointer histogramFilter = ImageToHistogramType::New() ;
      histogramFilter->SetInput( image ) ;
      histogramFilter->SetNumberOfBins( (unsigned int)minMaxFilter->GetMaximum() ) ;
      histogramFilter->Compute() ;
      typedef itk::OtsuMultipleThresholdsCalculator< typename ImageToHistogramType::HistogramType > ThresholdFilterType ;
      typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New() ;
      thresholdFilter->SetInputHistogram( histogramFilter->GetOutput() ) ;
      thresholdFilter->SetNumberOfThresholds( numberThreshold ) ;
      thresholdFilter->Update() ;
      typename ThresholdFilterType::OutputType thresholds = thresholdFilter->GetOutput() ;
      std::cout << "Computed threshold(s) (before scaling): " ;
      for( int i = 0 ; i < numberThreshold ; i++ )
      {
         std::cout << thresholds[ i ] << " " ;
      }
      std::cout << std::endl ;
      if( numberThreshold == 1 )
      {
         lowerThreshold = thresholds[ 0 ] ;
         upperThreshold = minMaxFilter->GetMaximum () ;
      }
      else
      {
         lowerThreshold = 1 ;
         upperThreshold = thresholds[ 0 ] ;
      }
    }
    if( threshold )
    {
       return EXIT_SUCCESS ;
    }
    ///////Grayscale erosion///////////////////////////////
    typedef itk::Neighborhood< T , 3 > NeighborhoodType ;
    NeighborhoodType neighborhood ;
    neighborhood.SetRadius( 2 ) ;
    for( int i = -2 ; i <= 2 ; i++ )
    {
       for( int j = -2 ; j <= 2 ; j++ )
       {
          for( int k = -2 ; k <= 2 ; k++ )
          {
             itk::Offset< 3 > offset ;
             offset[ 0 ] = i ;
             offset[ 1 ] = j ;
             offset[ 2 ] = k ;
             int val = ( i*i + j*j + k*k <= 2*2 ? 1 : 0 ) ;
             neighborhood[ offset ] = val ;
          }
       }
}
   /* typedef itk::FlatStructuringElement< 3 > StructuringElementType ;
    itk::Size< 3 > size ;
    size.Fill( 2 ) ;
    StructuringElementType structuringElement = StructuringElementType::Ball( size ) ;*/
    typedef itk::GrayscaleErodeImageFilter< InputImageType ,
                                            InputImageType ,
                                            NeighborhoodType
                                                  > GrayscaleErodeFilterType ;
    /*typedef itk::GrayscaleErodeImageFilter< InputImageType ,
                                            InputImageType ,
                                            StructuringElementType
                                                  > GrayscaleErodeFilterType ;*/
    for( int i = 0 ; i < 2 ; i++ )
    {
      typename GrayscaleErodeFilterType::Pointer grayscaleErode = GrayscaleErodeFilterType::New() ;
      grayscaleErode->SetKernel( neighborhood );
      //grayscaleErode->SetKernel( structuringElement );
      grayscaleErode->SetInput( image ) ;
      grayscaleErode->Update() ;
      image = grayscaleErode->GetOutput() ;
      image->DisconnectPipeline() ;
    }
    upperThreshold *= scaleThreshold ;
    lowerThreshold *= scaleThreshold ;
    if( upperThreshold > std::numeric_limits< T >::max() )
    {
       upperThreshold = std::numeric_limits< T >::max() ;
    }
    if( upperThreshold < std::numeric_limits< T >::min() )
    {
       upperThreshold = std::numeric_limits< T >::min() ;
    }
    if( lowerThreshold > std::numeric_limits< T >::max() )
    {
       lowerThreshold = std::numeric_limits< T >::max() ;
    }
    if( lowerThreshold < std::numeric_limits< T >::min() )
    {
       lowerThreshold = std::numeric_limits< T >::min() ;
    }
    //////Threshold///////////////////////////////////////////
    typedef itk::Image< unsigned char, 3 > MaskImageType ;
    typedef itk::BinaryThresholdImageFilter< InputImageType , MaskImageType > ThresholdFilterType ;
    typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New() ;
    thresholdFilter->SetUpperThreshold( (T)upperThreshold ) ;
    thresholdFilter->SetLowerThreshold( (T)lowerThreshold ) ;
    thresholdFilter->SetInput( image ) ;
    thresholdFilter->Update() ;
    ///////Fill holes///////////////////////////////////////////////////
    typedef itk::GrayscaleFillholeImageFilter< MaskImageType , MaskImageType > FillHoleFilterType ;
    FillHoleFilterType::Pointer fillHoleFilter = FillHoleFilterType::New() ;
    fillHoleFilter->SetInput( thresholdFilter->GetOutput() ) ;
    fillHoleFilter->Update() ;
    //////Binary erode////////////////////////////////////////
    MaskImageType::Pointer mask = fillHoleFilter->GetOutput() ;
    mask->DisconnectPipeline() ;
    typedef itk::BinaryErodeImageFilter< MaskImageType ,
                                         MaskImageType ,
                                         NeighborhoodType
    > ErodeFilterType ;
    /*typedef itk::BinaryErodeImageFilter< MaskImageType ,
                                         MaskImageType ,
                                         StructuringElementType
                                         > ErodeFilterType ;*/
    for( int i = 0 ; i < numberOfErosions ; i++ )//6
    {
      typename ErodeFilterType::Pointer binaryErode = ErodeFilterType::New() ;
      //binaryErode->SetKernel( structuringElement ) ;
      binaryErode->SetKernel( neighborhood ) ;
      binaryErode->SetErodeValue( 255 ) ;
      binaryErode->SetInput( mask ) ;
      binaryErode->Update() ;
      mask = binaryErode->GetOutput() ;
      mask->DisconnectPipeline() ;
    }
    /////////////////Connected components
    typedef itk::Image< unsigned int , 3 > CCImageType ;
    typedef itk::ConnectedComponentImageFilter< MaskImageType , CCImageType > CCFilterType ;
    CCFilterType::Pointer ccFilter = CCFilterType::New() ;
    ccFilter->SetInput( mask ) ;
    ccFilter->Update() ;
    /////////////////Remove everything except the largest region
    typedef itk::CastImageFilter< CCImageType , MaskImageType > CastFilterType ;
    CastFilterType::Pointer castFilter = CastFilterType::New() ;
    castFilter->SetInput( ccFilter->GetOutput() ) ;
    castFilter->Update() ;
    mask = castFilter->GetOutput() ;
    itk::ImageRegionIterator< MaskImageType > it( mask , mask->GetLargestPossibleRegion() ) ;
    std::vector< long > veccount ;
    for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
    {
      if( it.Get() >= veccount.size() )
      {
        veccount.push_back( 1 ) ;
      }
      else
      {
        veccount[ it.Get() ]++ ;
      }
    }
    long index = 0 ;
    long count = 0 ;
    for( size_t i = 1 ; i < veccount.size() ; i++ )
    {
      if( veccount[ i ] > count )
      {
        index = i ;
        count = veccount[ i ] ;
      }
    }
    for( it.GoToBegin() ; !it.IsAtEnd() ; ++it )
    {
      it.Set( it.Get() != index ? 0 : 1 ) ;
    }
    /////////////////Binary Dilate /////////////////////////
    typedef itk::BinaryDilateImageFilter< MaskImageType ,
                                          MaskImageType ,
                                          NeighborhoodType
                                                > DilateFilterType ;
   /* typedef itk::BinaryDilateImageFilter< MaskImageType ,
                                          MaskImageType ,
                                          StructuringElementType
                                                > DilateFilterType ;*/
    for( int i = 0 ; i < numberOfErosions + 2 ; i++ )//8
    {
      typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New() ;
      //binaryDilate->SetKernel( structuringElement ) ;
      binaryDilate->SetKernel( neighborhood ) ;
      binaryDilate->SetDilateValue( 1 ) ;
      binaryDilate->SetInput( mask ) ;
      binaryDilate->Update() ;
      mask = binaryDilate->GetOutput() ;
      mask->DisconnectPipeline() ;
    }
    ////////////////Fill holes //////////////////////////////
    fillHoleFilter->SetInput( mask ) ;
    fillHoleFilter->Update() ;
    fillHoleFilter->GetOutput() ;
    ///////////////////    Save output volume   ////////////////////////////////////////////////
    typedef itk::ImageFileWriter< MaskImageType > ImageWriterType ;
    typename ImageWriterType::Pointer writer = ImageWriterType::New() ;
    writer->SetFileName( outputVolume ) ;
    writer->SetInput( fillHoleFilter->GetOutput() ) ;
    writer->SetUseCompression( true ) ;
    writer->Update() ;
    return EXIT_SUCCESS ;
}

void GetImageType( std::string fileName,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType
                 )
{
  typedef itk::Image< unsigned char , 3 > ImageType ;
  itk::ImageFileReader< ImageType >::Pointer imageReader 
              = itk::ImageFileReader< ImageType >::New() ;
  imageReader->SetFileName( fileName.c_str() ) ;
  imageReader->UpdateOutputInformation() ;
  pixelType = imageReader->GetImageIO()->GetPixelType() ;
  componentType = imageReader->GetImageIO()->GetComponentType() ;
}



int TemplateInputVolume( std::string inputVolume , int argc , char * argv[] )
{
  itk::ImageIOBase::IOPixelType pixelType ;
  itk::ImageIOBase::IOComponentType componentType ;
  try
  {
    GetImageType ( inputVolume , pixelType , componentType ) ;
    // This filter handles all image component types
    switch( componentType )
    {
      case itk::ImageIOBase::UCHAR:
        return DoIt< unsigned char >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::CHAR:
        return DoIt< char >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::USHORT:
        return DoIt< unsigned short >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::SHORT:
        return DoIt< short >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::UINT:
        return DoIt< unsigned int >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::INT:
        return DoIt< int >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::ULONG:
        return DoIt< unsigned long >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::LONG:
        return DoIt< long >( argc , argv ) ;
        break ;
      case itk::ImageIOBase::FLOAT:
      case itk::ImageIOBase::DOUBLE:
        std::cerr << "Input image voxels must be integers " << std::endl ;
        return EXIT_FAILURE ;
        break ;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
        default:
        return EXIT_FAILURE ;
        break ;
    }
  }
  catch( itk::ExceptionObject &excep )
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl ;
    std::cerr << excep << std::endl ;
    return EXIT_FAILURE ;
  }
}


int main( int argc, char * argv[] )
{
  /////////////////   Get arguments using GenerateCLP
  //// parser to get the input volume filename  ///////////////////////
  PARSE_ARGS;
  if( !outputVolume.compare("") && !threshold )
  {
     std::cerr << "You must specify an output file name" << std::endl ;
     return EXIT_FAILURE ;
  }
  if( outputVolume.compare("") && threshold )
  {
     std::cout << "Warning: No file will be written out" << std::endl ;
  }
  /////////////////   Read input volume data type and instantiate
  //// the corresponding templated filter  function    ////////////////////
  return TemplateInputVolume( inputVolume , argc , argv );
}
