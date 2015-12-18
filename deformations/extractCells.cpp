// Standard library
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <fstream>
#include <cstdlib>

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
#include <DGtal/io/viewers/Viewer3D.h>
#include <DGtal/topology/ExplicitDigitalSurface.h>


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
  bool exportEvolver = false;

  // Program options
  namespace po = boost::program_options;
  po::options_description general_opt( "Allowed options are: " );
  general_opt.add_options()
    ("help,h",      "display this message")
    ("dimension,d", po::value< std::vector<unsigned int> >(), "dimensions of the image.")
    ("implicit,i",  po::value< std::string >(), "raw (double) real image where the level-set of value 0 represents the cell interfaces.")
    ("thickness,t", po::value< Real >( &thickness )->default_value( thickness ), "the thickening parameter for the implicit surface.")
    ("priority,p",  po::value< bool >( &calcPriority )->default_value( calcPriority ), "control if the priority is calculated (if a real image is given).")
    ("scale,s",     po::value< DGtal::uint32_t >(&priorityScale)->default_value(priorityScale), "Factor applied to the implicit data to get the cell priority.")
    ("label,l",     po::value< std::string >(), "raw (unsigned short int) label image. Used to add labels border and to better calculate priority.")
    ("view,v", po::value< std::string >()->default_value( "Normal" ), "specifies if the surface is viewed as is (Normal) or if places close to singularities are highlighted (Singular), or if unsure places should not be displayed (Hide), or if no view is wanted (no)." )
    ("evolver,e",   po::value< std::string >(), "if set, the result is exported to Surface Evolver format with the given file name." )
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
         KSpace::periodic
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
  // TODO
  if ( labelImage != nullptr )
    {
      trace.beginBlock( "Separating surface from labels." );

      for ( Point pt : domain )
        {
          const Cell spel = K.uSpel(pt);

          for ( Dimension d = 0; d < dimension; ++d )
            {
              const Cell nextSpel = K.uIncident( spel, d, true );

              if ( (*labelImage)(pt) != (*labelImage)( K.uCoords(nextSpel) ) )
                fullComplex.insertCell( nextSpel, unsureData );
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
  // Collapsing cells
  trace.beginBlock( "Collapsing cells." );
  functions::ccops::collapse(
      fullComplex,
      fullComplex.begin(), fullComplex.end(),
      typename CC::DefaultCellMapIteratorPriority{},
      true, true, true
  );

  trace.info() << "       K     = " << fullComplex << endl;
  trace.endBlock();
  trace.info() << endl;

  // Viewing the result.
  if ( vm[ "view" ].as<std::string>() != "no" )
    {

      //-------------- Copy to close Khalimsky space -------------------------------------------
      trace.beginBlock( "Copying to a closed Khalimsky space." );

      KSpace cK;
      cK.init(
              domain.lowerBound(),
              domain.upperBound(),
              KSpace::closed
      );
      CC fullClosedComplex(cK);

      for ( auto const& cell : fullComplex )
        fullClosedComplex.insertCell( cell );

      for ( DGtal::Dimension i = 0; i < KSpace::dimension; ++i )
        {
          for ( KSpace::Integer x = cK.lowerCell().myCoordinates[(i+1)%3]; x != cK.upperCell().myCoordinates[(i+1)%3]+1; ++x )
            {
              for ( KSpace::Integer y = cK.lowerCell().myCoordinates[(i+2)%3]; y != cK.upperCell().myCoordinates[(i+2)%3]+1; ++y )
                {
                  Cell p;
                  p.myCoordinates[(i+1)%3] = x;
                  p.myCoordinates[(i+2)%3] = y;

                  p.myCoordinates[i] = cK.lowerCell().myCoordinates[i];
                  if ( fullComplex.belongs( K.uCell( p.myCoordinates ) ) )
                    fullClosedComplex.insertCell( p );

                  p.myCoordinates[i] = cK.upperCell().myCoordinates[i];
                  if ( fullComplex.belongs( K.uCell( p.myCoordinates ) ) )
                    fullClosedComplex.insertCell( p );
                }
            }
        }
      trace.info() << "     C K     = " << fullClosedComplex << endl;

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
      for ( auto it = fullClosedComplex.begin( 0 ), itEnd = fullClosedComplex.end( 0 ); it != itEnd; ++it, ++idx )
        {
          Cell cell = it->first;
          indices[ cell ] = idx;
          points.push_back( cK.uKCoords(cell)/2 - Point::diagonal(1) );
          mesh.addVertex(   cK.uKCoords(cell)/2 - Point::diagonal(1) );
        }

      for ( auto it = fullClosedComplex.begin( 2 ), itEnd = fullClosedComplex.end( 2 ); it != itEnd; ++it )
        {
          Cell cell = it->first;
          bool fixed = it->second.data & CC::FIXED;

          Cells bdry = fullClosedComplex.cellBoundary( cell, true );
          std::vector<unsigned int> face_idx;
          for ( auto itC = bdry.begin(), itCE = bdry.end(); itC != itCE; ++itC )
            if ( fullClosedComplex.dim( *itC ) == 0 )
              face_idx.push_back( indices[*itC] );

          /*
          if ( ( ! fixed ) && hide ) continue;
          Color color = highlight
            ? ( fixed ? Color(128,255,128) : Color::White )
            : Color::White;
          */
          Color color = Color::White;
          if ( hide ) color.alpha(64);

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
        }
      trace.endBlock();

      //-------------- View surface -------------------------------------------
      const Point lowerKCoords = cK.uKCoords( cK.lowerCell() );
      const Point upperKCoords = cK.uKCoords( cK.upperCell() );

      QApplication application(argc,argv);
      Viewer3D<Space,KSpace> viewer( K );
      viewer.setWindowTitle("simple Volume Viewer");
      viewer.show();
      viewer << mesh;
      
      // Display lines that are not in the mesh.
      for ( auto it = fullClosedComplex.begin( 1 ), itE = fullClosedComplex.end( 1 ); it != itE; ++it )
        {
          Cell cell  = it->first;
          bool fixed  = it->second.data & CC::FIXED;
          std::vector<Cell> dummy;
          std::back_insert_iterator< std::vector<Cell> > outIt( dummy );
          fullClosedComplex.directCoFaces( outIt, cell );

          const auto coords = cK.uKCoords(cell);
          bool isOnBorder = false;
          for ( std::size_t i = 0; i < dimension & !isOnBorder; ++i )
            isOnBorder |= ( coords[i] == lowerKCoords[i] || coords[i] == upperKCoords[i] );

          //if ( ! dummy.empty() && !( isOnBorder && fixed ) )     continue;
          if ( ! dummy.empty() ) continue;

          Cells bdry = fullClosedComplex.cellBoundary( cell, true );
          Cell v0    = *(bdry.begin() );
          Cell v1    = *(bdry.begin() + 1);
          //if ( ( ! fixed ) && hide ) continue;
          Color color = highlight || hide ? Color::Red : Color::White;
          /*
          Color color = highlight
            ? ( fixed ? Color::White : Color(128,255,128) )
            : Color::White;
          */
          viewer.setLineColor( color );
          viewer.addLine( points[ indices[ v0 ] ], points[ indices[ v1 ] ], 2.0 );
        }
      
      // Display points that are not in the mesh.
      for ( auto it = fullClosedComplex.begin( 0 ), itE = fullClosedComplex.end( 0 ); it != itE; ++it )
        {
          Cell cell  = it->first;
          std::vector<Cell> dummy;
          std::back_insert_iterator< std::vector<Cell> > outIt( dummy );
          fullClosedComplex.directCoFaces( outIt, cell );

          if ( ! dummy.empty() ) continue;

          Color color = highlight || hide ? Color::Red : Color::White;
          
          viewer.setLineColor( color );
          viewer.addBall( points[ indices[ cell ] ], 2.0 );
        }

      viewer << Viewer3D<Space,KSpace>::updateDisplay;
      application.exec();
    }
  // End viewing.


  //-------------- View surface -------------------------------------------
  if ( vm.count( "evolver" ) )
    {
      const std::string fileName = vm[ "evolver" ].as<std::string>();

      trace.beginBlock( "Surface Evolver export to " + fileName );
      
      std::ofstream fileStream( fileName+".fe", std::ofstream::out | std::ofstream::binary );
      fileStream.imbue( std::locale() ); // Dot separator for decimal numbers
      fileStream  << "// extracCells.cpp\n";
      fileStream  << "TORUS_FILLED\n\n";
      fileStream  << "periods\n"
                  << "1 0 0\n" 
                  << "0 1 0\n" 
                  << "0 0 1\n\n";

      //-----------------------------------------------------------------------
      trace.beginBlock( "Indexing and writing vertices" );
      fileStream  << "vertices\n";
      std::map< Cell, std::size_t > index0;
      std::size_t idx = 1;
      for ( auto it = fullComplex.begin( 0 ), itEnd = fullComplex.end( 0 ); it != itEnd; ++it, ++idx )
        {
          index0[ it->first ] = idx;
          fileStream << idx;
          
          const Point coords = K.uKCoords( it->first ); 
          
          for ( Dimension i = 0; i < dimension; ++i )
            fileStream << " " << ( double(coords[i]) / (2*extent[i]+1) );

          fileStream << "\n";
        }
      fileStream << "\n";
      trace.endBlock();
      
      //-----------------------------------------------------------------------
      trace.beginBlock( "Indexing and writing edges" );
      fileStream  << "edges\n";
      std::map< Cell, std::size_t > index1;
      idx = 1;
      for ( auto it = fullComplex.begin( 1 ), itEnd = fullComplex.end( 1 ); it != itEnd; ++it, ++idx )
        {
          const Cell edge = it->first;
          index1[ edge ] = idx;
          fileStream << idx;

          const Dimension dir = *( K.uDirs( edge ) );
          const Cell vertex1  = K.uIncident( edge, dir, false );
          const Cell vertex2  = K.uIncident( edge, dir, true );
          const Point coords1 = K.uKCoords( vertex1 );
          const Point coords2 = K.uKCoords( vertex2 );

          fileStream  << " " << index0[ vertex1 ]
                      << " " << index0[ vertex2 ];
          
          for ( std::size_t i = 0; i < dimension; ++i )
            {
              if ( std::abs( coords1[i] - coords2[i] ) <= 2 )
                fileStream << " *";
              else
                fileStream << ( coords1[i] < coords2[i] ? " -" : " +" );
            }

          fileStream << "\n";
        }
      fileStream << "\n";
      index0.clear();
      trace.endBlock();
      
      //-----------------------------------------------------------------------
      trace.beginBlock( "Indexing faces" );

      // Stores index and usage of each sign
      struct FaceInfo
        {
          std::size_t idx;
          std::size_t side[2];
        };

      std::map< Cell, FaceInfo > index2;
      idx = 1;
      for ( auto it = fullComplex.begin( 2 ), itEnd = fullComplex.end( 2 ); it != itEnd; ++it, ++idx )
        {
          index2[ it->first ] = { idx, { 0, 0 } };
        }

      trace.endBlock();

      //-----------------------------------------------------------------------
      trace.beginBlock( "Indexing bodies and writing faces" );

      fileStream << "faces\n";

      // Storing bodies to count them before writing (needed for volume calculation).
      std::vector< std::vector< long int > > bodies;

      // Predicate for the faces.
      struct FacePredicate
        {
          using Surfel = KSpace::SCell;
          KSpace const* K;
          CC const* myComplex;
          inline bool operator() ( Surfel const& cell ) const
            {
              return myComplex->belongs( 2, K->unsigns(cell) );
            }
        };

      const FacePredicate predicate{ &K, &fullComplex };
      const SurfelAdjacency<dimension> adjacency( true );

      while ( ! index2.empty() )
        {
          std::vector< long int > faces;
          const auto startFace = index2.begin();
          const KSpace::SCell startCell = K.signs( startFace->first, startFace->second.side[0] == 0 ? K.NEG : K.POS ); // index2 doesn't contain any face with the two sides used.

          trace.info() << "Body #" << ( bodies.size() + 1 ) << " starting at: " << startCell << endl;

          const auto surface = DGtal::ExplicitDigitalSurface< KSpace, FacePredicate >( K, predicate, adjacency, startCell, true );
          for ( auto const& cell : surface )
            {
              auto it = index2.find( K.unsigns( cell ) );
              ASSERT( it != index2.end() );

              faces.push_back( it->second.idx * ( K.sSign(cell) == K.POS ? 1 : -1 ) );
              it->second.side[ K.sSign(cell) == K.POS ? 1 : 0 ] = bodies.size() + 1; // Tag this side as used.

              if ( it->second.side[ K.sSign(cell) == K.POS ? 0 : 1 ] != 0 ) // It the other side is also used, write and erase this face.
                {
                  const Cell face = it->first;
                  fileStream << it->second.idx;

                  const Dimension d = K.uOrthDir( face );
                  // TODO: why is the orientation for d=1 not what Surface Evolver expects ?
                  if ( d != 1 )
                    {
                      fileStream  << " -" << index1[ K.uIncident( face, (d+1)%3, true ) ]
                        << " -" << index1[ K.uIncident( face, (d+2)%3, false ) ]
                        << "  " << index1[ K.uIncident( face, (d+1)%3, false  ) ]
                        << "  " << index1[ K.uIncident( face, (d+2)%3, true  ) ]
                        << " frontcolor " << ( ( it->second.side[1] - 1 ) % 15 + 1 )
                        << " backcolor "  << ( ( it->second.side[0] - 1 ) % 15 + 1 )
                        << "\n";
                    }
                  else
                    {
                      fileStream  << "  " << index1[ K.uIncident( face, (d+1)%3, true ) ]
                        << " -" << index1[ K.uIncident( face, (d+2)%3, true  ) ]
                        << " -" << index1[ K.uIncident( face, (d+1)%3, false  ) ]
                        << "  " << index1[ K.uIncident( face, (d+2)%3, false ) ]
                        << " frontcolor " << ( ( it->second.side[1] - 1 ) % 15 + 1 )
                        << " backcolor "  << ( ( it->second.side[0] - 1 ) % 15 + 1 )
                        << "\n";
                    }
                  
                  index2.erase( it );
                }
            }

          bodies.emplace_back( std::move(faces) );
        }
      
      fileStream << "\n";
      index1.clear();
      trace.endBlock();
      
      //-----------------------------------------------------------------------
      trace.beginBlock( "Writing bodies" );

      fileStream << "bodies\n";
      for ( std::size_t i = 0; i < bodies.size(); ++i )
        {
          fileStream << (i+1);
          for ( const auto idx : bodies[i] )
            fileStream << " " << idx;
          fileStream << " volume 1/" << bodies.size() << "\n";
        }

      trace.endBlock();
      
      //-----------------------------------------------------------------------
      fileStream << 
        "read\n"
        "hessian_normal\n"
        "gogo := { g 5; V; r; g 5; r; g 5; convert_to_quantities; hessian; hessian; }\n";

      fileStream.close();
      trace.endBlock();
    }


  return 0;
}
