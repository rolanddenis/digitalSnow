#include <sstream>
#include <iomanip>
#include <cstddef>
#include <map>
#include <numeric>
#include <cmath>
#include <limits>
#include <random>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

/////////////////////
#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>

// Dimension (default 3)
#ifndef DIMENSION
  #define DIMENSION 3
#endif
// See: http://stackoverflow.com/questions/1489932/c-preprocessor-and-token-concatenation
#define PASTER(dim) Z ## dim ## i
#define EVALUATOR(dim) PASTER(dim)
#define ZNI EVALUATOR( DIMENSION )

using namespace DGtal;
using namespace ZNI;
using namespace std;

// Useful functions
#include "FunctorConstImage.h"

// Massive Multi Phase-Field for dripping
#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"
#include "MultiPhaseFieldDrip.h"

// IO functions
#include <DGtal/io/writers/VTKLightWriter.h>
#include <DGtal/io/writers/RawWriter.h>
#if   DIMENSION == 2
  #include "deformationDisplay2d.h"
  #include "DGtal/io/readers/PGMReader.h"
#elif DIMENSION == 3
  #include "DGtal/io/readers/VolReader.h"
#endif

// Dimension as a variable
namespace
{
  static const
  unsigned int dimension = DIMENSION ;
}

///////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{

  DGtal::trace.info() << dimension << "d dripping using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;

  // Default options
  size_t dsize = 64;          // Domain size
  double tstep = 0.25;        // Time step
  size_t disp_step  = 1;      // Display step
  size_t max_step = 1;        // Maximum number of steps
  double epsilon = 3.;        // Interface width
  size_t max_phase_cnt = 64;  // Maximal number of phases
  unsigned int seed = std::random_device{}(); // Seed for the random number generator.

  string outputFiles  = "interface";  // Output files basename

#if   DIMENSION == 2
    string outputFormat = "raster"; // Output files format
#elif DIMENSION == 3
    string outputFormat = "vol";    // Output files format
#endif

  // Command-line options description
  namespace po = boost::program_options;
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h",          "display this message")
    ("domainSize,d",    po::value<size_t>(&dsize)->default_value(dsize), "Domain size (if default starting interface)" )
    ("maxPhaseCnt,p",   po::value<size_t>(&max_phase_cnt)->default_value(max_phase_cnt),
        "Maximal number of phases." )
    ("timeStep,t",      po::value<double>(&tstep)->default_value(tstep), "Time step for the evolution" )
    ("displayStep",     po::value<size_t>(&disp_step)->default_value(disp_step), "Number of time steps between 2 drawings" )
    ("stepsNumber,n",   po::value<size_t>(&max_step)->default_value(max_step), "Maximal number of steps" )
    ("epsilon,e",       po::value<double>(&epsilon)->default_value(epsilon), "Interface width (only for phase fields)" )
    ("seed",            po::value<unsigned int>(&seed), "Seed used to initialize the random number generator.")
    ("outputFiles,o",   po::value<string>(&outputFiles)->default_value(outputFiles), "Output files basename" )
#if   DIMENSION == 2
    ("outputFormat,f", po::value<string>(&outputFormat)->default_value(outputFormat),
     "Output files format: either <raster> (image, default) or <vector> (domain representation)" );
#elif DIMENSION == 3
    ("outputFormat,f",  po::value<string>(&outputFormat)->default_value(outputFormat),
        "Output files format: either <png> (3d to 2d with QGLViewer), <pngc> (3d to 2d with Cairo) or <vol> (3d)" );
#endif

  // Command-line parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);
  po::notify(vm);
  if(vm.count("help") || argc<=1)
    {
      trace.info()<< "Evolution of a " << dimension << "d dripping" << std::endl
        << "Basic usage: "<<std::endl
        << argv[0] << " [other options] -t <time step> -n <number of steps>"
        << std::endl
        << general_opt << "\n";
      return 1;
    }

  /////////////////////////////////////////////////////////////////////////////
  // Options validity
  // Files format
#if   DIMENSION == 2
  if ( (outputFormat != "vector") && (outputFormat != "raster") )
    {
      trace.info() << "format is expected to be either <vector> or <raster> " << std::endl;
      return 1;
    }
#elif DIMENSION == 3
  if ( (outputFormat != "png") && (outputFormat != "vol") && (outputFormat != "pngc") )
    {
      trace.info() << "format is expected to be either <png>, <pngc> or <vol> " << std::endl;
      return 1;
    }
#endif

  // Checking epsilon
  if (epsilon <= 0)
  {
      trace.error() << "epsilon should be greater than 0" << std::endl;
      return 1;
  }

  /////////////////////////////////////////////////////////////////////////////

  // Random number generator
  trace.info() << "Seed used for the random number generator: " << seed << std::endl;
  std::mt19937 gen(seed);

  // Domain
  const Point p = Point::diagonal(0);
  const Point q = Point::diagonal(dsize-1); // Domain size = q-p+1
  const Domain domain = Domain( p, q );
  trace.info() << "Domain = " << domain << std::endl << std::endl;

  // Initial label image.
  using Label = short unsigned int;
  using LabelImage =  ImageContainerBySTLVector<Domain, Label>;
  LabelImage labelImage( domain );


  // Field image
  typedef ImageContainerBySTLVector<Domain, double> FieldImage;

  // ApproximatedMultiImage
  using real = double;

  using LabelledMap = DGtal::LabelledMap<real, 64, long unsigned int, 4, 4 >;
  //using LabelledMap = DGtal::BigLabelledMap<real, (1ul<<7)-1, 2, 10>;

  using Approximation = DGtal::approximations::NegativeTolValueApproximation<real>;
  //using Approximation = DGtal::approximations::NoValueApproximation<real>;
  //using BoundingBox = AxisAlignedBoundingBox< Domain, unsigned int>;
  using BoundingBox = NoBoundingBox< Domain >;
  using ApproximatedMultiImage = DGtal::ApproximatedMultiImage<Domain, LabelledMap, Approximation, BoundingBox>;

  // Multi phase-field
  epsilon /= labelImage.extent()[0];
  tstep = epsilon*epsilon;

  MultiPhaseFieldDrip< FieldImage, ApproximatedMultiImage >
      evolver( domain, max_phase_cnt, epsilon );

  evolver.updateLabels( labelImage );

  DGtal::trace.beginBlock( "Dripping" );

  // Initial state export
  std::stringstream s;
  s << outputFiles << setfill('0') << std::setw(6) << 0;
#if   DIMENSION == 2
  drawContours( labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
  //writePartition( labelImage, s.str(), outputFormat );
#endif
  
  auto const implicitImage = makeFunctorConstImage( labelImage.domain(),
      [&evolver, epsilon] ( Point const& aPoint ) -> real
      {
        real max1 = 0, max2 = 0;
        for ( auto value : evolver.getPhasesContainer()(aPoint) )
        {
          if ( value.second >= max1 ) { max2 = max1; max1 = value.second; }
          else if ( value.second > max2 ) { max2 = value.second; }
        }
        //return max1 - max2;
        return 2 * epsilon * std::atanh( std::min(max1 - max2, static_cast<real>(1 - 1e-8) ) );
      }
  );

  auto const storageImage = makeFunctorConstImage ( labelImage.domain(),
      [&evolver] ( Point const& aPoint )
      {
        return static_cast<unsigned int>(evolver.getPhasesContainer()(aPoint).size());
      }
  );

  // VTK export
    {
      VTKLightWriter<Domain> vtk(s.str(), labelImage.domain(),
          evolver.myRealExtent / ( labelImage.domain().upperBound() - labelImage.domain().lowerBound() + Point::diagonal() ) );
      /*
         for (size_t j = 0; j < evolver.getNumPhase(); ++j)
         {
         stringstream s_phase;
         s_phase << "phi" << setfill('0') << std::setw(2) << j;
         vtk << s_phase.str() << evolver.getPhase(j);
         }
         */
      vtk << "label"    << labelImage;
      vtk << "implicit" << implicitImage;
      vtk << "storage"  << storageImage;
    }

  // Informations
  evolver.dispInfos();

  // Time integration
  double sumt = 0;
  for (unsigned int i = 1; i <= max_step; ++i)
    {
      DGtal::trace.info() << "iteration # " << i << std::endl;

      // Update
      trace.beginBlock("Iteration");
      evolver.update( tstep );
      trace.endBlock();

      // Display
      if ( (i % disp_step) == 0 )
        {
          // Update labels
          const std::size_t label_cnt = evolver.updateLabels( labelImage );

          // Display phase field informations
          evolver.dispInfos();

          // Export
          trace.beginBlock("Export");



          std::stringstream s;
          s << outputFiles << setfill('0') << std::setw(6) << i;
#if   DIMENSION == 2
          drawContours( labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
          //writePartition( labelImage, s.str(), outputFormat );
#endif

          RawWriter< decltype(implicitImage) >::exportRaw<real>( s.str()+".imp.raw", implicitImage );
          RawWriter< LabelImage >::exportRaw<unsigned short int>( s.str()+".lab.raw", labelImage );

          // VTK export
          VTKLightWriter<Domain> vtk(s.str(), labelImage.domain(),
            evolver.myRealExtent / ( labelImage.domain().upperBound() - labelImage.domain().lowerBound() + Point::diagonal() ) );
          /*
             for (size_t j = 0; j < evolver.getNumPhase(); ++j)
             {
             stringstream s_phase;
             s_phase << "phi" << setfill('0') << std::setw(2) << j;
             vtk << s_phase.str() << evolver.getPhase(j);
             }
             */
          vtk << "label"    << labelImage;
          vtk << "implicit" << implicitImage;
          vtk << "storage"  << storageImage;
          vtk.close();

          // Volume of each phase
          /*
             DGtal::trace.info() << ( dimension == 2 ? "Area: " : "Volume: " );
             for (size_t j = 0; j < evolver.getNumPhase(); ++j)
             {
             DGtal::trace.info() << "V(" << j << ") = " << getVolume<double>(evolver.getPhase(j));
             }
             DGtal::trace.info() << std::endl;
             */

          trace.endBlock();

          /*
          evolver.updateDomainSize();

          std::cout << "myRealExtent = " << evolver.myRealExtent << std::endl;
          for ( Dimension j = 0; j < dimension; ++j )
            {
              const auto re = evolver.myRealExtent[j];
              if ( std::isnan( re ) || re <= 0.2 || re >= 5.0 )
                {
                  std::cerr << "Error: invalid domain size !!!" << std::endl;
                  return 1;
                }
            }
          */

          if ( label_cnt <= 0.00001 * disp_step * labelImage.domain().size() * ( evolver.getNumPhase() == max_phase_cnt ? 0.1 : 1. ) )
            if ( ! evolver.addPhase( gen ) )
              break;

        }

      sumt += tstep;

      // Volume of each phase
      /*
         typedef std::map<typename LabelImage::Value, unsigned int> Histo;
         Histo histo;
         calcHistogram( *labelImage, histo );
         DGtal::trace.info() << "Volume: ";
         for ( Histo::const_iterator it = histo.begin(); it != histo.end(); ++it )
         DGtal::trace.info() << "V(" << it->first << ") = " << it->second;
         DGtal::trace.info() << std::endl;
         */


      DGtal::trace.info() << "Time spent: " << sumt << std::endl << std::endl;
    }

  DGtal::trace.endBlock();

  if ( dimension == 3)
  {
      std::cout << "Command line to extract cells:" << std::endl;
      std::cout << "extractCells" << std::setprecision(20)
          << " -d " << labelImage.extent()[0]
          << " -d " << labelImage.extent()[1]
          << " -d " << labelImage.extent()[2]
          << " -S " << evolver.myRealExtent[0]
          << " -S " << evolver.myRealExtent[1]
          << " -S " << evolver.myRealExtent[2]
          //<< " -l " << ( s.str() + ".lab.raw" )
          //<< " -v no"
          << std::endl;
  }

  return 0;
}


/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
