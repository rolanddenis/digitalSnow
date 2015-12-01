// Standard library
#include <iostream>
#include <string>
#include <map>

// Boost
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

// DGtal
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/kernel/SpaceND.h>
#include <DGtal/kernel/domains/HyperRectDomain.h>
#include <DGtal/topology/KhalimskySpaceND.h>
#include <DGtal/topology/CubicalComplex.h>
#include <DGtal/topology/CubicalComplexFunctions.h>
#include <DGtal/images/ImageContainerBySTLVector.h>
#include <DGtal/io/readers/GenericReader.h>

#define DIMENSION 3

using namespace DGtal;

int main ( int argc, char* argv[] )
{
  // Aliases
  using Real    = double;
  using Integer = int;
  const Integer dimension = DIMENSION;
  using Space   = SpaceND< dimension, Integer >;
  using Point   = Space::Point;
  using Domain  = HyperRectDomain< Space >;
  
  using KSpace  = KhalimskySpaceND< dimension, Integer >;
  using Cell    = KSpace::Cell;
  using CCMap   = std::map< Cell, CubicalCellData >;
  using CC      = CubicalComplex< KSpace, CCMap >;

  using Label   = unsigned int;
  using RealImage   = ImageContainerBySTLVector<Domain, Real>;
  using LabelImage  = ImageContainerBySTLVector<Domain, Label>;
  
  // Other aliases
  using std::cout;
  using std::cerr;
  using std::endl;

  // Default parameters
  bool periodic = true;
  Real tickness = 0;

  // Program options
  namespace po = boost::program_options;
  po::options_description general_opt( "Allowed options are: " ).add_options()
    ("help,h",      "display this message")
    ("periodic,p",  po::value< bool >( &periodic )->default_value( periodic ), "define if the domain is periodic.")
    ("implicit,i",  po::value< std::string >(), "real image where the level-set of value 0 represents the cell interfaces.")
    ("thickness,t", po::value< Real >( &thickness )->default_value( thickness ), "the thickening parameter for the implicit surface.")
  ;

  // Parsing program options
  bool isParseOK = true;
  po::variables_map vm;
  try
    {
      po::store( po::parse_command_line( argc, argv, general_opt ), vm );
    }
  catch ( const std::exception& ex )
    {
      isParseOK = false;
      cerr << "Error checking program options: " << ex.what() << endl;
    }
  po::notify( vm );

  // Displaying help
  if ( ! IsParseOK || vm.count("help") || ! vm.count("implicit") )
    {
      cerr  << "Usage: " << argv[0] << " -i <file> [options]" << endl
            << general_opt << endl;
      return 1;
    }

  // Default space and domain.
  KSpace K;
  Domain domain;

  // Reading real image
  RealImage realImage( Domain() );
  if ( vm.count("implicit") )
    {
      trace.beginBlock( "Reading real image." );
      const std::string realImageName = vm["implicit"].as<std::string>();
      realImage = GenericReader< Domain >::import( realImageName );
      domain = realImage.domain();
      trace.endBlock();
    }

  // Initializing Khalimsky space
  if ( periodic )
    {
      K.init( 
             domain.lowerBound() - Point::diagonal(1),
             domain.upperBound() + Point::diagonal(1),
             true
      );
    }
  else
    {
      K.init( domain.lowerBound(), domain.upperBound(), true );
    }
  
  // Initializing cellular complex.
  CC fullComplex( K );
  CubicalCellData unsureData( 0 );
  CubicalCellData sureData( CC::FIXED );

  // Using interface implicit representation
  if ( vm.count("implicit") )
    {
      trace.beginBlock( "Filling cellular complex with thickened interface from implicit representation." );
      
      for ( Point pt : domain )
        {
          Cell spel = K.uSpel(pt);
          if ( realImage(pt) <= thickness )
            fullComplex.insertCell( spel, unsureData );
        }
      trace.endBlock();

      trace.info() << "    K = " << fullComplex << endl;
      fullComplex.close();
      trace.info() << " C1 K = " << fullComplex << endl;

    }

  return 0;
}
