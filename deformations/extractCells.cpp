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
#include <boost/function_output_iterator.hpp>

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
#include "DGtal/io/viewers/Viewer3D.h"

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
  using RealPoint = Space::RealPoint;

  using KSpace  = KhalimskySpaceND< dimension, Integer >;
  using Cell    = KSpace::Cell;
  using Cells   = KSpace::Cells;
  using CCMap   = std::map< Cell, CubicalCellData >;
  using CC      = CubicalComplex< KSpace, CCMap >;
  using CellMapIterator = CC::CellMapIterator;
  using CCData  = CC::Data;

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
    ("label,l",     po::value< std::string >(), "raw (unsigned short int) label image. Used to add labels border and better calculate priority.")
    ("view,v", po::value< std::string >()->default_value( "Normal" ), "specifies if the surface is viewed as is (Normal) or if places close to singularities are highlighted (Singular), or if unsure places should not be displayed (Hide)." )
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

  /////////////////////////////////////////////////////////////////////////////
  // Reading real image
  RealImage* realImage = nullptr;
  if ( vm.count("implicit") )
    {
      trace.beginBlock( "Reading real image." );
      const std::string realImageName = vm["implicit"].as<std::string>();
      realImage = new RealImage( RawReader< RealImage >::importRaw<Real>( realImageName, extent ) );
      trace.endBlock();
      trace.info() << endl;
    }

  /////////////////////////////////////////////////////////////////////////////
  // Reading label image
  LabelImage* labelImage = nullptr;
  if ( vm.count("label") )
    {
      trace.beginBlock( "Reading label image." );
      const std::string labelImageName = vm["label"].as<std::string>();
      labelImage = new LabelImage( RawReader< LabelImage >::importRaw<Label>( labelImageName, extent ) );
      trace.endBlock();
      trace.info() << endl;
    }

  // Initializing Khalimsky space
  K.init(
         domain.lowerBound(),
         domain.upperBound(),
         true
  );

  // Initializing cellular complex.
  CC fullComplex( K );
  CubicalCellData unsureData( 0 );
  CubicalCellData sureData( CC::FIXED );

  /////////////////////////////////////////////////////////////////////////////
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

      trace.info() << "           K = " << fullComplex << endl;
      fullComplex.close();
      trace.info() << "        C1 K = " << fullComplex << endl;

      trace.endBlock();
      trace.info() << endl;
    }

  /////////////////////////////////////////////////////////////////////////////
  // Using labels to add missed surface.
  if ( labelImage != nullptr )
    {
      trace.beginBlock( "Separating surface from labels." );

      for ( Point pt : domain )
        {
          const Cell spel = K.uSpel(pt);

          for ( Dimension d = 0; d < dimension; ++d )
            {
              Point nextPt = pt + Point::base(d);
              Point periodicNextPt = nextPt;
              if ( periodicNextPt[d] == extent[d] )
                periodicNextPt[d] = 0;
              if ( (*labelImage)(pt) != (*labelImage)(periodicNextPt) )
                fullComplex.insertCell( K.uIncident( spel, d, true ), unsureData );
            }
        }

      trace.info() << " C1 K +    S = " << fullComplex << endl;
      fullComplex.close();
      trace.info() << " C1 K + C1 S = " << fullComplex << endl;

      trace.endBlock();
      trace.info() << endl;
    }

  /////////////////////////////////////////////////////////////////////////////
  // Computes priority
  if ( calcPriority && realImage != nullptr )
    {
      ASSERT( dimension > 0 );

      trace.beginBlock( "Computing priority" );

      for ( std::size_t i = 0; i <= dimension; ++i )
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
      trace.info() << endl;
    }

  /////////////////////////////////////////////////////////////////////////////
  // Fixing boundary cells
  trace.beginBlock( "Fixing boundary cells." );
  const Point lowerKCoords = K.uKCoords( K.lowerCell() );
  const Point upperKCoords = K.uKCoords( K.upperCell() );
  cout << lowerKCoords << " ; " << upperKCoords << endl;

  std::size_t bdryCnt = 0, innerCnt = 0;
  functions::ccops::filterCellsWithinBounds (
      fullComplex,
      lowerKCoords, upperKCoords,
      boost::make_function_output_iterator( // Boundary cells
          [&fullComplex, &bdryCnt] ( Cell const& cell ) { fullComplex.findCell( cell )->second.data |= CC::FIXED; ++bdryCnt; }
      ),
      boost::make_function_output_iterator( [&innerCnt] ( Cell ) { ++innerCnt; } ) // Inner cells
  );
  trace.info() << innerCnt << " inner cells and " << bdryCnt << " boundary cells." << endl;
  trace.endBlock();
  trace.info() << endl;

  /////////////////////////////////////////////////////////////////////////////
  // Collapsing inner cells
  trace.beginBlock( "Collapsing inner cells." );

  functions::ccops::collapse(
      fullComplex,
      fullComplex.begin(), fullComplex.end(),
      typename CC::DefaultCellMapIteratorPriority{},
      true, true, true
  );

  trace.info() << "       K     = " << fullComplex << endl;
  trace.endBlock();
  trace.info() << endl;

  /////////////////////////////////////////////////////////////////////////////
  // Fixing inner cells and freeing boundary cells.
  trace.beginBlock( "Fixing inner cells and freeing boundary cells." );
  
  innerCnt = 0;
  std::vector<Cell> bdryCells;

  functions::ccops::filterCellsWithinBounds (
      fullComplex,
      lowerKCoords, upperKCoords,
      boost::make_function_output_iterator( // Boundary cells
          [&fullComplex, &bdryCells] ( Cell const& cell ) { fullComplex.findCell(cell)->second.data &= ~CC::FIXED; bdryCells.push_back(cell); }
      ),
      boost::make_function_output_iterator( // Inner cells
          [&fullComplex, &innerCnt] ( Cell const& cell ) { fullComplex.findCell( cell )->second.data |= CC::FIXED; ++innerCnt; }
      )
  );
  trace.info() << innerCnt << " inner cells and " << bdryCells.size() << " boundary cells." << endl;
  trace.endBlock();
  trace.info() << endl;

  /////////////////////////////////////////////////////////////////////////////
  // Fixing boundary cells incident to inner cells.
  trace.beginBlock( "Fixing boundary cells incident to inner cells." );
  bdryCnt = 0;
  for ( auto const& cell : bdryCells )
    {
      if ( K.uDim(cell) == dimension-2 )
        {
          const Point ptCoords = K.uKCoords(cell);
          for ( std::size_t i = 0; i < dimension; ++i )
            {
              if ( ptCoords[i] == lowerKCoords[i] || ptCoords[i] == upperKCoords[i] )
                {
                  const auto incCell = fullComplex.findCell( dimension-1, K.uIncident( cell, i, ptCoords[i] == lowerKCoords[i] ) );
                  if ( incCell != fullComplex.end(dimension-1) && ( incCell->second.data & CC::FIXED ) )
                    {
                      fullComplex.findCell( dimension-2, cell )->second.data |= CC::FIXED;
                      ++bdryCnt;
                      break;
                    }
                }
            }
        }
    }
  trace.info() << bdryCnt << " boundary cells fixed." << endl;
  trace.endBlock();
  trace.info() << endl;

  /////////////////////////////////////////////////////////////////////////////
  // Duplicating boundary cells by periodicity
  trace.beginBlock( "Duplicating boundary cells by periodicity." );
  std::vector< std::pair<Cell,CCData> > toBeInserted;

  for ( auto const& cell : bdryCells )
    {
      std::vector<Point> cellClones = { K.uKCoords(cell) };
      // Generates clones coordinates.
      for ( std::size_t d = 0; d < dimension; ++d )
        {
          if ( cellClones[0][d] == lowerKCoords[d] || cellClones[0][d] == upperKCoords[d] )
            {
              const auto coordShift = ( cellClones[0][d] == lowerKCoords[d] ? 2 : -2 ) * extent[d];
              for ( std::size_t i = 0, size = cellClones.size(); i < size; ++i )
                cellClones.push_back( cellClones[i] + Point::base( d, coordShift ) );
            }
        }

      // Generates clones cell & data.
      const auto it = fullComplex.findCell( cell );
      for ( std::size_t i = 1; i < cellClones.size(); ++i )
        toBeInserted.emplace_back( Cell(cellClones[i]), it->second );
    }

  // Insert clones
  trace.info() << toBeInserted.size() << " boundary cells clones by periodicity." << endl;
  for ( auto const& element : toBeInserted )
    {
      const auto it = fullComplex.findCell( element.first );
      if ( it != fullComplex.end( K.uDim( element.first ) ) )
        it->second.data |= element.second.data & CC::FIXED;
      else
        fullComplex.insertCell( element.first, element.second );
    }

  toBeInserted.clear();

  trace.info() << "       K + P = " << fullComplex << endl;
  trace.endBlock();
  trace.info() << endl;

  /////////////////////////////////////////////////////////////////////////////
  // Collapsing boundary cells
  trace.beginBlock( "Collapsing boundary cells." );
  
  bdryCells.clear();
  innerCnt = 0;
  
  functions::ccops::filterCellsWithinBounds (
      fullComplex,
      lowerKCoords, upperKCoords,
      std::back_inserter( bdryCells ), // Boundary cells
      boost::make_function_output_iterator( [&innerCnt] ( Cell ) { ++innerCnt; } ) // Inner cells
  );
  trace.info() << innerCnt << " inner cells and " << bdryCells.size() << " boundary cells." << endl;
  
  functions::ccops::collapse(
      fullComplex,
      bdryCells.cbegin(), bdryCells.cend(),
      typename CC::DefaultCellMapIteratorPriority{},
      true, true, true
  );

  trace.info() << "       K     = " << fullComplex << endl;
  trace.endBlock();
  trace.info() << endl;
  

  //-------------- Create Mesh -------------------------------------------
  trace.beginBlock( "Create Mesh. " );
  std::string view = vm[ "view" ].as<std::string>();
  bool highlight = ( view == "Singular" );
  bool hide      = ( view == "Hide" );
  Mesh<Point> mesh( true );
  std::map<Cell, unsigned int> indices;
  std::vector<Point> points;
  int idx = 0;
  for ( auto it = fullComplex.begin( 0 ), itEnd = fullComplex.end( 0 ); it != itEnd; ++it, ++idx )
    {
      Cell cell = it->first;
      indices[ cell ] = idx;
      points.push_back( K.uKCoords(cell)/2 - Point::diagonal(1) );
      mesh.addVertex(   K.uKCoords(cell)/2 - Point::diagonal(1) );
    }

  for ( auto it = fullComplex.begin( 2 ), itEnd = fullComplex.end( 2 ); it != itEnd; ++it )
    {
      Cell cell = it->first;
      bool fixed = it->second.data & CC::FIXED;

      Cells bdry = fullComplex.cellBoundary( cell, true );
      std::vector<unsigned int> face_idx;
      for ( auto itC = bdry.begin(), itCE = bdry.end(); itC != itCE; ++itC )
        {
          if ( fullComplex.dim( *itC ) == 0 )
            face_idx.push_back( indices[ *itC ] );
        }

      if ( ( ! fixed ) && hide ) continue;
      Color color = highlight
        ? ( fixed ? Color(128,255,128) : Color::White )
        : Color::White;
      Vector diag03 = points[ face_idx[ 0 ] ] - points[ face_idx[ 3 ] ];
      Vector diag12 = points[ face_idx[ 1 ] ] - points[ face_idx[ 2 ] ];
      if ( diag03.dot( diag03 ) <= diag12.dot( diag12 ) )
        {
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 3 ], color );
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 3 ], face_idx[ 2 ], color );
        }
      else
        {
          mesh.addTriangularFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 2 ], color );
          mesh.addTriangularFace( face_idx[ 1 ], face_idx[ 3 ], face_idx[ 2 ], color );
        }
      //mesh.addQuadFace( face_idx[ 0 ], face_idx[ 1 ], face_idx[ 3 ], face_idx[ 2 ], color );
     }
  trace.endBlock();

  //-------------- View surface -------------------------------------------
  QApplication application(argc,argv);
  Viewer3D<Space,KSpace> viewer( K );
  viewer.setWindowTitle("simple Volume Viewer");
  viewer.show();
  viewer << mesh;
  // Display lines that are not in the mesh.
  for ( auto it = fullComplex.begin( 1 ), itE = fullComplex.end( 1 ); it != itE; ++it )
    {
      Cell cell  = it->first;
      bool fixed  = it->second.data & CC::FIXED;
      std::vector<Cell> dummy;
      std::back_insert_iterator< std::vector<Cell> > outIt( dummy );
      fullComplex.directCoFaces( outIt, cell );

      const auto coords = K.uKCoords(cell);
      bool isOnBorder = false;
      for ( std::size_t i = 0; i < dimension & !isOnBorder; ++i )
        isOnBorder |= ( coords[i] == lowerKCoords[i] || coords[i] == upperKCoords[i] );

      if ( ! dummy.empty() && !( isOnBorder && fixed ) )     continue;

      Cells bdry = fullComplex.cellBoundary( cell, true );
      Cell v0    = *(bdry.begin() );
      Cell v1    = *(bdry.begin() + 1);
      if ( ( ! fixed ) && hide ) continue;
      Color color = highlight
        ? ( fixed ? Color::White : Color(128,255,128) )
        : Color::White;
      viewer.setLineColor( color );
      viewer.addLine( points[ indices[ v0 ] ], points[ indices[ v1 ] ], 1/2.0 );
    }
  viewer << Viewer3D<Space,KSpace>::updateDisplay;
  return application.exec();
}
