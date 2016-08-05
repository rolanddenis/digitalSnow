#include <sstream>
#include <iomanip>
#include <cstddef>
#include <map>
#include <numeric>
#include <cmath>
#include <limits>

#include <vector>
#include <set>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;

/////////////////////
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

#define DIMENSION 3

using namespace DGtal;
using namespace Z3i;
using namespace std;

// Evolvers

// Massive Multi Phase-Field
#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"
#include "MultiPhaseField2.h"

// Useful functions
#include "deformationFunctions.h"

// IO functions
#include "VTKWriter.h"
#include <DGtal/io/writers/RawWriter.h>
#include <DGtal/io/writers/VolWriter.h>
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/RawReader.h"

// Dimension as a variable
namespace
{
  static const constexpr
  unsigned int dimension = DIMENSION ;
}

///////////////////////////////////////////////////////////////////////////////
template <
  typename TImage,
  typename TFormatList,
  typename TEvolver
>
void writePartition( TImage const & anImage, TEvolver const & anEvolver, std::string const & aFileName, TFormatList const & aFormatList )
{
  for ( auto const & fileFormat : aFormatList )
    {
      if      ( fileFormat == "vol" )
        VolWriter< TImage >::exportVol( aFileName + ".vol", anImage );
      else if ( fileFormat == "raw" )
        RawWriter< TImage >::exportRaw8( aFileName + ".raw", anImage );
      else if ( fileFormat == "vtk" )
        {
          VTKWriter< typename TImage::Domain > vtk( aFileName, anImage.domain() );
          vtk << "label" << anImage;
          for (size_t j = 0; j < anEvolver.getNumPhase(); ++j)
            {
              stringstream s_phase;
              s_phase << "phi" << setfill('0') << std::setw(2) << j;
              vtk << s_phase.str() << anEvolver.getPhase(j);
            }

        }
    }
}

template <
  typename TImage,
  typename TFormatList
>
void writePartition( TImage const & anImage, std::string const & aFileName, TFormatList const & aFormatList )
{
  for ( auto const & fileFormat : aFormatList )
    {
      if      ( fileFormat == "vol" )
        VolWriter< TImage >::exportVol( aFileName + ".vol", anImage );
      else if ( fileFormat == "raw" )
        RawWriter< TImage >::exportRaw8( aFileName + ".raw", anImage );
      else if ( fileFormat == "vtk" )
        {
          VTKWriter< typename TImage::Domain > vtk( aFileName, anImage.domain() );
          vtk << "label" << anImage;
        }
    }
}

///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{

  DGtal::trace.info() << dimension << "d snow interface evolution using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;

  // Default options
  double timeStep = 0.25;         // Time step
  size_t displayStep  = 1;        // Display step
  size_t maxStep = 1;             // Maximum number of steps
  bool   doNotCalcDist = false;   // Don't calculate the distance function to initiliaze phase field (for <massiveMultiPhaseField>)
  double epsilon = 3.;            // Interface width
  std::string outputFileName  = "interface";  // Output files basename
  std::string outputFormatString = "vol"; // Output files format.

  // Command-line options description
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h",          "display this message")
    ("domainSize,d",    po::value<string>(), "Domain size (if initializing with <raw> file format) in LxHxP format." )
    ("inputImage,i",    po::value<string>(), "Image to initialize the starting interface (<vol> or <raw> format)" )
    ("timeStep,t",      po::value<double>( &timeStep )->default_value( timeStep ), "Time step for the evolution. Shouldn't be greather than epsilon²." )
    ("displayStep",     po::value<size_t>( &displayStep )->default_value( displayStep ), "Number of time steps between 2 exports." )
    ("stepsNumber,n",   po::value<size_t>( &maxStep )->default_value( maxStep ), "Maximal number of steps." )
    ("noDist",          po::bool_switch( &doNotCalcDist ), "Do not initialize the phase field with distance function.")
    ("epsilon,e",       po::value<double>( &epsilon )->default_value( epsilon ), "Interface width as a multiple of space step." )
    ("outputFiles,o",   po::value<string>( &outputFileName )->default_value( outputFileName ), "Output files basename" )
    ("outputFormat,f",  po::value<string>( &outputFormatString )->default_value( outputFormatString ),
        "Output files format list (comma separeted) composed from either <vol>, <raw> or <vtk>." );

  // Command-line parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);
  po::notify(vm);

  if( vm.count("help") || argc <= 1 )
    {
      trace.info()<< "Evolution of a " << dimension << "d snow interface" << std::endl
        << "Basic usage: "<<std::endl
        << argv[0] << " [other options] -t <time step> -n <number of steps>"
        << std::endl
        << general_opt << "\n";
      return 1;
    }

  // Parsing files format.
  std::set< std::string > outputFormat;

  {
    std::vector< std::string > outputFormatList;
    boost::algorithm::split(  outputFormatList,
                              vm["outputFormat"].as< std::string >(),
                              boost::algorithm::is_any_of( "," ) );

    outputFormat.insert( outputFormatList.cbegin(), outputFormatList.cend() );
  }

  // Checking output files format.
  if ( outputFormat.size() == 0 )
    {
      trace.error() << "You must choose at least one output file format." << std::endl;
      return 1;
    }

  for ( auto const & fileFormat : outputFormat )
    {
      if ( (fileFormat != "vol") && (fileFormat != "vtk") && (fileFormat != "raw") && (fileFormat != "") )
        {
          trace.error() << "Output format is expected to be either <vol>, <raw> or <vtk>." << std::endl;
          return 1;
        }
    }

  // Loading image.
  if ( vm.count("inputImage") == 0 )
    {
      trace.error() << "You must specify an input image file." << std::endl;
      return 1;
    }

  using Label = unsigned char;
  using LabelImage = ImageContainerBySTLVector<Domain, Label>;
  LabelImage labelImage( Domain{} );

  {
    const auto imageFileName = vm["inputImage"].as<std::string>();

    // Getting file extension.
    const auto fileExtension = boost::filesystem::path( imageFileName ).extension();

    // Reading file depending on his extension.
    if ( fileExtension == ".vol" )
      {
        trace.beginBlock( "Reading image \"" + imageFileName + "\" in vol format..." );
        labelImage = VolReader<LabelImage>::importVol( imageFileName );
        trace.endBlock();
      }
    else if ( fileExtension == ".raw" || fileExtension == "" )
      {
        trace.beginBlock( "Reading image \"" + imageFileName + "\" in raw format..." );

        if ( vm.count("domainSize") == 0 )
          {
            trace.error() << "When importing a raw file, you must specify the image size." << std::endl;
            return 1;
          }

        // Parsing domain size.
        Vector domainSize;
        std::size_t dimCounter = 0;

        std::vector< std::string > domainSizeString;
        boost::algorithm::split(  domainSizeString,
                                  vm["domainSize"].as< std::string >(),
                                  boost::algorithm::is_any_of( "x" ) );

        if ( domainSizeString.size() > domainSize.size() )
          {
            trace.error() << "Too many dimensions specified for image size." << std::endl;
            return 1;
          }

        // Converting each token from string to integer.
        for ( const auto & token : domainSizeString )
          try
            {
              domainSize[ dimCounter++ ] = std::stoi( token );
            }
          catch( std::invalid_argument const & exception )
            {
              trace.error() << "Invalid domain size " << token << std::endl;
              return 1;
            }

        // Completing non specified dimensions.
        for ( ; dimCounter < domainSize.size(); ++dimCounter )
          domainSize[ dimCounter ] = domainSize[ dimCounter-1 ];

        // Displaying domain size.
        trace.info() << "Domain size is " << domainSize[0] << "x"
                                          << domainSize[1] << "x"
                                          << domainSize[2] << "."
                                          << std::endl;

        // Reading file.
        labelImage = DGtal::RawReader<LabelImage>::importRaw8( imageFileName, domainSize );

        trace.endBlock();
      }
    else
      {
        trace.error() << "Unknow file extension " << fileExtension << " for reading \"" << imageFileName << "\"." << std::endl;
        return 1;
      }
  }

  // Epsilon's validity
  if (epsilon <= 0)
    {
      trace.error() << "epsilon should be greater than 0" << std::endl;
      return 1;
    }

  // Domain
  const Domain domain = labelImage.domain();
  trace.info() << std::endl << "Domain = " << domain << std::endl << std::endl;

  // Field image
  typedef ImageContainerBySTLVector<Domain, double> FieldImage;

  // ApproximatedMultiImage
  using real = double;
  using LabelledMap = DGtal::LabelledMap<real, 64, long unsigned int, 1, 2 >;
  using Approximation = DGtal::approximations::NegativeTolValueApproximation<real>;
  //using BoundingBox = AxisAlignedBoundingBox< Domain, unsigned int>;
  using BoundingBox = NoBoundingBox< Domain >;
  using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation, BoundingBox>;

  MultiPhaseField2< LabelImage, FieldImage, ApproximatedMultiImage > evolver(labelImage, epsilon, !doNotCalcDist);

  trace.info() << std::endl;
  trace.beginBlock( "Deformation (massive multi phase field)" );

  // Initial state export
  evolver.updateLabels();
  std::stringstream s;
  s << outputFileName << setfill('0') << std::setw(4) << 0;
  writePartition( labelImage, evolver, s.str(), outputFormat );

  // Informations
  evolver.dispInfos();
  std::cout << std::endl;

  // Time integration
  double sumt = 0;
  for (unsigned int i = 1; i <= maxStep; ++i)
    {
      DGtal::trace.info() << "iteration # " << i << std::endl;

      // Update
      trace.beginBlock("Iteration");
      sumt += evolver.update( timeStep );
      trace.endBlock();
      trace.info() << std::endl;

      // Display
      if ( (i % displayStep) == 0 )
        {
          // Update labels
          const std::size_t label_cnt = evolver.updateLabels();

          // Display phase field informations
          evolver.dispInfos();
          std::cout << std::endl;

          // Export
          trace.beginBlock("Export");

          std::stringstream s;
          s << outputFileName << setfill('0') << std::setw(4) << (i/displayStep);
          writePartition( labelImage, evolver, s.str(), outputFormat );

          trace.endBlock();

        }

      DGtal::trace.info() << "Time spent: " << sumt << std::endl << std::endl;
    }

  DGtal::trace.endBlock();

  return 0;

}

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */

