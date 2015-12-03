// Standard library
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>

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
#include <DGtal/io/readers/RawReader.h>

#define DIMENSION 3

using namespace DGtal;

template <
  typename TKhalimsky,
  typename TRealImage,
  typename TLabelImage
>
std::pair< typename TRealImage::Value, typename TLabelImage::Value >
interpValue ( TKhalimsky const& K, TRealImage const* realImage, TLabelImage const* labelImage, typename TKhalimsky::Cell const& cell )
{
  using Value   = typename TRealImage::Value;
  using Label   = typename TLabelImage::Value;
  using Domain  = typename TRealImage::Domain;
  using Dimension = typename Domain::Dimension;
  using Point   = typename Domain::Point;

  assert( realImage != nullptr );

  if ( K.uDim(cell) == Domain::dimension ) // Get values from images.
    {
      const Domain domain = realImage->domain();
      Point pt = ( K.uKCoords(cell) - Point::diagonal(1) ) / 2;
      for ( std::size_t i = 0; i < Domain::dimension; ++i )
        {
          if ( pt[i] == domain.lowerBound()[i]-1 )
            pt[i] = domain.upperBound()[i];
          else if ( pt[i] == domain.upperBound()[i]+1 )
            pt[i] = domain.lowerBound()[i];
        }

      if ( labelImage != nullptr )
        return { (*realImage)(pt), (*labelImage)(pt) };
      else
        return { (*realImage)(pt), 0 };
    }
  else if ( K.uDim(cell) == Domain::dimension-1 ) // Mean of neighbor taking into account the labels.
    {
      const Dimension d = K.uOrthDir(cell);
      auto interp1 = interpValue( K, realImage, labelImage, K.uIncident( cell, d, false ) );
      auto interp2 = interpValue( K, realImage, labelImage, K.uIncident( cell, d, true ) );
      if ( interp1.second == interp2.second )
        {
          interp1.first = ( interp1.first + interp2.first ) / 2;
        }
      else
        {
          interp1.first = ( interp1.first - interp2.first ) / 2;
          if ( interp1.first < 0 )  interp1.second = interp2.second;
          interp1.first = std::abs( interp1.first );
        }

      return interp1;
    }
  else // Mean of neighbor values.
    {
      std::size_t cnt = 0;
      typename TRealImage::Value value = 0;
      for ( auto const& incCell : K.uUpperIncident( cell ) )
        {
          value += interpValue( K, realImage, labelImage, incCell ).first;
          ++cnt;
        }

      return { value/cnt, 0 };
    }
}


int main ( int argc, char* argv[] )
{
  // Aliases
  using Real    = double;
  using Integer = int;
  const Integer dimension = DIMENSION;
  using Space   = SpaceND< dimension, Integer >;
  using Point   = Space::Point;
  using Domain  = HyperRectDomain< Space >;
  using Dimension = Space::Dimension;
  using Vector  = Space::Vector;
  
  using KSpace  = KhalimskySpaceND< dimension, Integer >;
  using Cell    = KSpace::Cell;
  using CCMap   = std::map< Cell, CubicalCellData >;
  using CC      = CubicalComplex< KSpace, CCMap >;
  using CellMapIterator = CC::CellMapIterator;

  using Label   = unsigned short int;
  using RealImage   = ImageContainerBySTLVector<Domain, Real>;
  using LabelImage  = ImageContainerBySTLVector<Domain, Label>;
  
  // Other aliases
  using std::cout;
  using std::cerr;
  using std::endl;

  // Default parameters
  Real thickness = 0;
  bool calcPriority = true;
  DGtal::uint32_t priorityScale = 100000;


  // Program options
  namespace po = boost::program_options;
  po::options_description general_opt( "Allowed options are: " );
  general_opt.add_options()
    ("help,h",      "display this message")
    ("dimension,d", po::value< std::vector<unsigned int> >(), "dimensions of the image.")
    ("implicit,i",  po::value< std::string >(), "raw (double) real image where the level-set of value 0 represents the cell interfaces.")
    ("thickness,t", po::value< Real >( &thickness )->default_value( thickness ), "the thickening parameter for the implicit surface.")
    ("priority,p",  po::value<bool>( &calcPriority )->default_value( calcPriority ), "control if the priority is calculated (if a real image is given).")
    ("scale,s",     po::value< DGtal::uint32_t >(&priorityScale)->default_value(priorityScale), "Factor applied to the implicit data to get the cell priority.")
    ("label,l",     po::value< std::string >(), "raw (unsigned short int) label image.")
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

  // Verifying dimensions
  if ( ! vm.count("dimension") || vm["dimension"].as< std::vector<unsigned int> >().size() != dimension )
    {
      cerr << dimension << " dimensions of the image must be specified." << endl;
      isParseOK = false;
    }

  // Displaying help
  if ( ! isParseOK || vm.count("help") || ( ! vm.count("implicit") && ! vm.count("label") ) )
    {
      cerr  << "Usage: " << argv[0] << " -i/-l <file> [options]" << endl
            << general_opt << endl;
      return 1;
    }

  // Domain
  Vector extent;
  const auto dimensions = vm["dimension"].as< std::vector<unsigned int> >();
  for ( std::size_t i = 0; i < dimension; ++i )
    extent[i] = dimensions[i];
  Domain domain ( Point::diagonal(0), extent - Point::diagonal(1) );
  trace.info() << "Domain = " << domain << endl;

  // Default space and domain.
  KSpace K;

  // Reading real image
  RealImage* realImage = nullptr;
  if ( vm.count("implicit") )
    {
      trace.beginBlock( "Reading real image." );
      const std::string realImageName = vm["implicit"].as<std::string>();
      realImage = new RealImage( RawReader< RealImage >::importRaw<Real>( realImageName, extent ) );
      trace.endBlock();
    }
  
  // Reading label image
  LabelImage* labelImage = nullptr;
  if ( vm.count("label") )
    {
      trace.beginBlock( "Reading label image." );
      const std::string labelImageName = vm["label"].as<std::string>();
      labelImage = new LabelImage( RawReader< LabelImage >::importRaw<Label>( labelImageName, extent ) );
      trace.endBlock();
    }

  // Initializing Khalimsky space
  K.init( 
         domain.lowerBound() - Point::diagonal(1),
         domain.upperBound() + Point::diagonal(1),
         true
  );
  
  // Initializing cellular complex.
  CC fullComplex( K );
  CubicalCellData unsureData( 0 );
  CubicalCellData sureData( CC::FIXED );

  // Using interface implicit representation
  if ( realImage != nullptr )
    {
      trace.beginBlock( "Filling cellular complex with thickened interface from implicit representation." );
      
      for ( Point pt : domain )
        {
          Cell spel = K.uSpel(pt);
          if ( (*realImage)(pt) <= thickness )
            fullComplex.insertCell( spel, unsureData );
        }

      trace.info() << "    K = " << fullComplex << endl;
      fullComplex.close();
      trace.info() << " C1 K = " << fullComplex << endl;
      
      trace.endBlock();
    }

  // Using labels to add missed surface.
  if ( labelImage != nullptr )
    {
      trace.beginBlock( "Separating surface from labels." );
      
      for ( Point pt : domain )
        {
          for ( Dimension d = 0; d < dimension; ++d )
            {
              Point nextPt = pt + Point::base(d);
              Point periodicNextPt = nextPt;
              if ( periodicNextPt[d] == extent[d] )   periodicNextPt[d] = 0;
              if ( (*labelImage)(pt) != (*labelImage)(periodicNextPt) )
                {
                  fullComplex.insertCell( Cell( Point::diagonal(1) + 2*pt + Point::base(d) ), unsureData );
                  // Periodicity (not sure it is needed)
                  //if ( periodicNextPt[d] == 0 )
                  //  fullComplex.insertCell( Cell( Point::diagonal(1) + 2*periodicNextPt - Point::base(d) ), unsureData );
                }
            }
        }
      
      trace.info() << " C1 K +    S = " << fullComplex << endl;
      fullComplex.close();
      trace.info() << " C1 K + C1 S = " << fullComplex << endl;

      trace.endBlock();
    }

  // Compute priority
  if ( calcPriority && realImage != nullptr )
    {
      ASSERT( dimension > 0 );
      
      trace.beginBlock( "Computing priority" );

      for ( std::size_t i = 0; i < dimension; ++i )
        {
          for ( auto it = fullComplex.begin(i), itEnd = fullComplex.end(i); it != itEnd; ++it )
            {
              const Cell cell = it->first;
              const auto value = static_cast<DGtal::uint32_t>( std::round( interpValue(K, realImage, labelImage, cell ).first * priorityScale ) );

              it->second.data &= ~CC::VALUE;
              it->second.data |= value & CC::VALUE;
            }
        }

      trace.endBlock();
    }

  // Getting boundary and inner cells
  trace.beginBlock( "Get boundary and inner cells." );
  std::vector<Cell> innerCells;
  std::vector<Cell> bndryCells;
  functions::ccops::filterCellsWithinBounds (
      fullComplex, 
      K.uKCoords( K3.lowerCell() ) + Point::diagonal(2), K3.uKCoords( K3.upperCell() ) - Point::diagonal(2),
      std::back_inserter( bdryCells ), std::back_inserter( innerCells ) 
  );
  trace.info() << innerCells.size() << " inner cells and " << bndryCells.size() << " boundary cells." << endl;
  trace.endBlock();
  
  // Collapsing inner cells


  return 0;
}
