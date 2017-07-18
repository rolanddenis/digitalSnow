#include <sstream>
#include <iomanip>
#include <cstddef>
#include <map>
#include <numeric>
#include <cmath>
#include <limits>

/////////////////////
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

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

// Evolvers
// Level-Set
#include "WeickertKuhneEvolver.h"

// Phase-Field
#include "ExactDiffusionEvolver.h"
#include "ExactReactionEvolver.h"
#include "ExplicitReactionEvolver.h"
#include "LieSplittingEvolver.h"
#include "MultiPhaseField.h"

// Local Level-Set
#include "SimplePointHelper.h"
#include "PartitionEvolver.h"

// Massive Multi Phase-Field
#include "AxisAlignedBoundingBox.h"
#include "ValueApproximations.h"
#include "ApproximatedMultiImage.h"
#include "MultiPhaseField2.h"


// Useful functions
#include "deformationFunctions.h"
#include "FunctorConstImage.h"

// IO functions
#include <DGtal/io/writers/VTKLightWriter.h>
#include <DGtal/io/writers/RawWriter.h>
#if   DIMENSION == 2
  #include "deformationDisplay2d.h"
  #include "DGtal/io/readers/PGMReader.h"
#elif DIMENSION == 3
  #include "deformationDisplay3d.h"
  #include "DGtal/io/readers/VolReader.h"
  #ifdef WITH_VISU3D_QGLVIEWER
    #include <QApplication> // Qt
  #endif
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

  DGtal::trace.info() << dimension << "d interface evolution using DGtal ";
  DGtal::trace.emphase() << "(version "<< DGTAL_VERSION << ")"<< std::endl;

#if DIMENSION == 3
  #ifdef WITH_VISU3D_QGLVIEWER
    // QApplication initialization with command-line parameters
    QApplication application(argc, argv);
  #endif
#endif

  // Default options
  size_t dsize = 64;          // Domain size
  double tstep = 0.25;        // Time step
  size_t disp_step  = 1;      // Display step
  size_t max_step = 1;        // Maximum number of steps
  string shape = "ball";      // Generated shape
  size_t phase_cnt = 2;       // Number of phases in each dimension (for mballs and rand initizaliation)
  string algo  = "massiveMultiPhaseField";  // Algorithm used for evolution
  bool   no_dist = false;     // Don't calculate the distance function to initiliaze phase field (for <massiveMultiPhaseField>)
  double balloon = 0.;        // Balloon Force
  double epsilon = 3.;        // Interface width
  bool flagWithCstVol = false;  // Volume Conservation
  size_t subSamplingRatio = 1;  // Sub-Sampling ratio
  size_t overSamplingRatio = 1; // Over-Sampling ratio
  size_t subRandomizeCount = 1; // New random labels within current labels.

  string outputFiles  = "interface";  // Output files basename

#if   DIMENSION == 2
    string outputFormat = "raster"; // Output files format
#elif DIMENSION == 3
    string outputFormat = "vol";    // Output files format
    bool flagWithVisu   = false;    // Interactive 3d visualization after evolution
#endif

  // Command-line options description
  po::options_description general_opt("Allowed options are");
  general_opt.add_options()
    ("help,h",          "display this message")
    ("inputImage,i",    po::value<string>(), "Binary image to initialize the starting interface (vol format)" )
    ("domainSize,d",    po::value<size_t>(&dsize)->default_value(dsize), "Domain size (if default starting interface)" )
    ("shape,s",         po::value<string>(&shape)->default_value(shape),
        "Generated shape: either <ball>, <flower>, <mballs> or <rand>" )
    ("seed",            po::value<unsigned int>(), "Seed used to initialize the random number generator for the rand shape.")
    ("phaseCnt,p",      po::value<size_t>(&phase_cnt)->default_value(phase_cnt),
        "Number of phases, for <mballs> and <rand> shapes." )
    ("timeStep,t",      po::value<double>(&tstep)->default_value(tstep), "Time step for the evolution" )
    ("displayStep",     po::value<size_t>(&disp_step)->default_value(disp_step), "Number of time steps between 2 drawings" )
    ("stepsNumber,n",   po::value<size_t>(&max_step)->default_value(max_step), "Maximal number of steps" )
    ("algo,a",          po::value<string>(&algo)->default_value(algo),
        "can be: \n <levelSet>  \n or <phaseField> \n or <multiPhaseField> \n or <massiveMultiPhaseField> \n or <localLevelSet>" )
    ("noDist",          po::bool_switch(&no_dist), "don't initialize the phase field with distance function, for <massiveMultiPhaseField>.")
    ("balloonForce,k",  po::value<double>(&balloon)->default_value(balloon), "Balloon force" )
    ("epsilon,e",       po::value<double>(&epsilon)->default_value(epsilon), "Interface width (only for phase fields)" )
    ("withCstVol",      po::bool_switch(&flagWithCstVol), "with volume conservation (only for phase fields)" )
    ("subSample",       po::value<size_t>(&subSamplingRatio)->default_value(subSamplingRatio), "sub-sampling ratio. Applied after image initialization but before over-sampling.")
    ("overSample",      po::value<size_t>(&overSamplingRatio)->default_value(overSamplingRatio), "over-sampling ratio. Applied after image initialization and over-sampling.")
    ("subRandomize",    po::value<size_t>(&subRandomizeCount)->default_value(subRandomizeCount), "Replaces each current label by the specified amount of new randomized label. Takes effect after image sampling.")
    ("outputFiles,o",   po::value<string>(&outputFiles)->default_value(outputFiles), "Output files basename" )
#if   DIMENSION == 2
    ("outputFormat,f", po::value<string>(&outputFormat)->default_value(outputFormat),
     "Output files format: either <raster> (image, default) or <vector> (domain representation)" );
#elif DIMENSION == 3
    ("outputFormat,f",  po::value<string>(&outputFormat)->default_value(outputFormat),
        "Output files format: either <png> (3d to 2d with QGLViewer), <pngc> (3d to 2d with Cairo) or <vol> (3d)" )
    ("withVisu",        po::bool_switch(&flagWithVisu), "Enables interactive 3d visualization after evolution" );
#endif

  // Command-line parsing
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, general_opt), vm);
  po::notify(vm);
  if(vm.count("help") || argc<=1)
    {
      trace.info()<< "Evolution of a " << dimension << "d interface" << std::endl
        << "Basic usage: "<<std::endl
        << argv[0] << " [other options] -t <time step> -n <number of steps>"
        << std::endl
        << general_opt << "\n";
      return 1;
    }

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

  // Generated shape
  if ( vm.count("inputImage") == 0 && shape != "ball" && shape != "flower" && shape != "mballs" && shape != "rand")
    {
      trace.info() << "if no input file is specified, shape is expected to be either <ball>, <flower>, <mballs> or <rand>" << std::endl;
      return 1;
    }

  // Evolution algorithm
  if ( algo != "levelSet" && algo != "phaseField" && algo != "multiPhaseField" && algo!= "massiveMultiPhaseField" && algo != "localLevelSet" )
    {
      trace.info() << "algo is expected to be either <levelSet>, <phaseField>, <multiPhaseField> or <localLevelSet> " << std::endl;
      return 1;
    }


  // Image and implicit function
  using Label = short unsigned int;
  typedef ImageContainerBySTLVector<Domain, Label> LabelImage;
  LabelImage* labelImage = NULL;

  if (vm.count("inputImage"))
    {
      string imageFileName = vm["inputImage"].as<std::string>();
      trace.emphase() << imageFileName <<std::endl;
      DGtal::trace.beginBlock("image reading...");

#if   DIMENSION == 2
      LabelImage tmp = PGMReader<LabelImage>::importPGM( imageFileName );
#elif DIMENSION == 3
      LabelImage tmp = VolReader<LabelImage>::importVol( imageFileName );
#endif

      labelImage = new LabelImage( tmp );
      DGtal::trace.endBlock();
    }
  else
    {
      DGtal::trace.beginBlock("image generating...");
      Point p = Point::diagonal(0);
      Point q = Point::diagonal(dsize-1); // Domain size = q-p+1
      Point c = Point::diagonal(dsize/2);
      labelImage = new LabelImage( Domain(p,q) );

      if ( (vm["shape"].as<std::string>()) == "flower" )
        initWithFlowerPredicate( *labelImage, c, (dsize*3/5)/2, (dsize*1/5)/2, 5 );
      else if ( (vm["shape"].as<std::string>()) == "ball" )
        initWithBallPredicate( *labelImage, c, (dsize*3/5)/2 );
      else if ( (vm["shape"].as<std::string>()) == "mballs" )
        initWithMultipleBalls( *labelImage, phase_cnt,  (dsize/phase_cnt) * 0.5, 1 );
      else
        {
          unsigned int seed;

          if ( vm.count("seed") )
            seed = initRandomly( *labelImage, phase_cnt, vm["seed"].as<unsigned int>() );
          else
            seed = initRandomly( *labelImage, phase_cnt );

          trace.info() << "Seed used to initialize image: " << seed << std::endl;
        }
      

      trace.info() << "starting interface initialized with a " << shape << std::endl;

      DGtal::trace.endBlock();
    }

  // Domain
  Domain d = Domain( labelImage->domain().lowerBound(), labelImage->domain().upperBound() );
  trace.info() << "Domain = " << d << std::endl << std::endl;

  // Sub-sampling
  if ( subSamplingRatio > 1 )
    {
      trace.beginBlock("Sub-sampling");

      const Point lower = d.lowerBound();
      const Point upper = lower + ( d.upperBound() - lower )/subSamplingRatio; // So that the domain becomes never empty.
      const Domain subDomain{ lower, upper };

      LabelImage* subLabelImage = new LabelImage{ subDomain };
      for ( auto point : subDomain )
        {
          auto const origPoint = lower + ( point - lower ) * subSamplingRatio;
          subLabelImage->setValue(point, (*labelImage)(origPoint));
        }

      d = subDomain;
      std::swap( labelImage, subLabelImage );
      delete subLabelImage;

      trace.info() << "Sub-sampled domain = " << d << std::endl;
      trace.endBlock();
      trace.info() << std::endl;
    }

  // Over-sampling
  if ( overSamplingRatio > 1 )
    {
      trace.beginBlock("Over-sampling");

      const Point lower = d.lowerBound();
      const Point upper = lower + ( d.upperBound() - lower + Point::diagonal(1) ) * overSamplingRatio - Point::diagonal(1);
      const Domain overDomain{ lower, upper };

      LabelImage* overLabelImage = new LabelImage{ overDomain };
      for ( auto point : overDomain )
        {
          auto const origPoint = lower + ( point - lower ) / overSamplingRatio;
          overLabelImage->setValue(point, (*labelImage)(origPoint));
        }

      d = overDomain;
      std::swap( labelImage, overLabelImage );
      delete overLabelImage;

      trace.info() << "Over-sampled domain = " << d << std::endl;
      trace.endBlock();
      trace.info() << std::endl;
    }

  // Add randomized labels within current labels.
  if ( subRandomizeCount > 1 )
    {
      trace.beginBlock( "Replacing current labels with new randomized labels." );

      trace.beginBlock( "Getting current labels" );

      std::set<Label> labels;
      for ( auto const& label : *labelImage )
        labels.insert( label );

      trace.endBlock();

      trace.beginBlock( "Adding new labels" );

      for ( auto label : labels )
        {
          trace.info() << "Label #" << label << std::endl;
          initRandomlyWithinLabel( *labelImage, label, subRandomizeCount );
        }

      trace.endBlock();

      trace.endBlock();
    }

#if DIMENSION == 3
  // Pre-visualization for choosing camera position
  if (outputFormat == "png")
    displayPartition( *labelImage );
#endif

  // Algorithm dispatch
  if ( algo == "levelSet" )
    {

      // Distance function
      ImageContainerBySTLVector<Domain,double> implicitFunction( d );
      initWithDT( *labelImage, implicitFunction );

#if DIMENSION == 2
      if (vm.count("withFunction"))
	      drawFunction( implicitFunction, vm["withFunction"].as<string>() );

      std::stringstream ss;
      ss << outputFiles << "0001";
      drawContour(implicitFunction, ss.str(), outputFormat);
#endif

      // Data functions
      ImageContainerBySTLVector<Domain,double> a( d );
      std::fill(a.begin(),a.end(), 1.0 );
      ImageContainerBySTLVector<Domain,double> b( d );
      std::fill(b.begin(),b.end(), 1.0 );
      ImageContainerBySTLVector<Domain,double> g( d );
      std::fill(g.begin(),g.end(), 1.0 );

      // Interface evolver
      WeickertKuhneEvolver<ImageContainerBySTLVector<Domain,double> > e(a, b, g, balloon, 1);

      DGtal::trace.beginBlock( "Deformation (Weickert's level set method)" );

      // Time integration
      double sumt = 0;
      for (unsigned int i = 1; i <= max_step; ++i)
        {
          DGtal::trace.info() << "iteration # " << i << std::endl;

          // Update
          e.update(implicitFunction, tstep);

          // Display
          if ((i % disp_step) == 0)
            {
              std::stringstream s;
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);
#if   DIMENSION == 2
	            drawContour(implicitFunction, s.str(), outputFormat);
#elif DIMENSION == 3
              updateLabelImage( *labelImage, implicitFunction );
              writePartition( *labelImage, s.str(), outputFormat );
#endif
            }

          sumt += tstep;

          // DGtal::trace.info() << "Area: " << getSize( implicitFunction ) << std::endl;
          //DGtal::trace.info() << "Volume: " << getSize(*labelImage, 0) << std::endl;

          DGtal::trace.info()
            << ( dimension == 2 ? "Area: " : "Volume: " )
            << getVolume<double>( implicitFunction )
            << std::endl;

          DGtal::trace.info() << "Time spent: " << sumt << std::endl;
        }

      //updateLabelImage( *labelImage, implicitFunction, 0 );
      //DGtal::trace.info() << "Volume: " << getSize(*labelImage, 0) << std::endl;
      DGtal::trace.endBlock();

#if DIMENSION == 3
      // Interactive display after the evolution
      if ( flagWithVisu )
        displayImageWithInfo( *labelImage, implicitFunction, a, b );
#endif

    }
  else if ( algo == "phaseField" )
    {
      // Option's validity
      if (epsilon <= 0)
        {
          trace.error() << "epsilon should be greater than 0" << std::endl;
          return 1;
        }

      // Distance function
      ImageContainerBySTLVector<Domain,double> implicitFunction( d );
      initWithDT( *labelImage, implicitFunction );

      // Computing the profile from the signed distance
      Profile p(epsilon);
      std::transform(implicitFunction.begin(), implicitFunction.end(), implicitFunction.begin(), p);

#if DIMENSION == 2
      if (vm.count("withFunction"))
        drawFunction( implicitFunction, vm["withFunction"].as<string>() );
#endif

      // Diffusion evolver
      typedef ExactDiffusionEvolver<ImageContainerBySTLVector<Domain,double> > Diffusion;
      Diffusion diffusion;

      // Exact reaction evolver (no possible volume conservation)
      /*
      typedef ExactReactionEvolver<ImageContainerBySTLVector<Domain,double> > Reaction;
      Reaction reaction( epsilon );
      */

      // Explicit reaction evolver (with possible volume conservation)
      typedef ExplicitReactionEvolver<
          ImageContainerBySTLVector<Domain,double>,
          ImageContainerBySTLVector<Domain,double>
        > Reaction;
      ImageContainerBySTLVector<Domain,double> a( Domain( implicitFunction.domain() ) );
      std::fill(a.begin(), a.end(), 0.0 );
      Reaction reaction( epsilon, a, balloon, flagWithCstVol );

      // Lie splitting
      LieSplittingEvolver<Diffusion,Reaction> e(diffusion, reaction);

      DGtal::trace.beginBlock( "Deformation (phase field)" );

      // Initial state export
      std::stringstream s;
      s << outputFiles << setfill('0') << std::setw(4) << 0;
#if DIMENSION == 2
      drawContour(implicitFunction, s.str(), outputFormat, 0.5);
#elif DIMENSION == 3
      writePartition( *labelImage, s.str(), outputFormat );
#endif

      // VTK export
        {
          VTKLightWriter<Domain> vtk(s.str(), implicitFunction.domain());
          vtk << "phi" << implicitFunction;
        }


      // Time integration
      double sumt = 0;
      for (unsigned int i = 1; i <= max_step; ++i)
        {
          DGtal::trace.info() << "iteration # " << i << std::endl;

          // Update
          e.update(implicitFunction, tstep);

          // Display
          if ((i % disp_step)==0)
            {
              std::stringstream s;
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);

#if   DIMENSION == 2
              drawContour(implicitFunction, s.str(), outputFormat, 0.5);
#elif DIMENSION == 3
              updateLabelImage( *labelImage, implicitFunction, 0.5 );
              writePartition( *labelImage, s.str(), outputFormat );
#endif

              // VTK export
              VTKLightWriter<Domain> vtk(s.str(), implicitFunction.domain());
              vtk << "phi" << implicitFunction;
            }

          sumt += tstep;

          DGtal::trace.info()
            << (dimension==2 ? "Area: " : "Volume: " )
            << getVolume<double>( implicitFunction )
            << std::endl;

          DGtal::trace.info() << "Time spent: " << sumt << std::endl;
        }

      // Post-processing
      // updateLabelImage( *labelImage, implicitFunction, 0.5 );
      // DGtal::trace.info() << "Volume: " << getSize(*labelImage, 0.5) << std::endl;

      DGtal::trace.endBlock();

#if DIMENSION == 3
      // Interactive display after the evolution
      if ( flagWithVisu )
        displayPartition( *labelImage );
#endif

    }
  else if ( algo == "multiPhaseField" )
    {
      // Options's validity
      if (epsilon <= 0)
        {
          trace.error() << "epsilon should be greater than 0" << std::endl;
          return 1;
        }

      // Field image
      typedef ImageContainerBySTLVector<Domain, double> FieldImage;

      // Diffusion evolver
      typedef ExactDiffusionEvolver< FieldImage > Diffusion;
      Diffusion diffusion;

      // Exact Reaction Evolver (no possible volume conservation)
      /*
      typedef ExactReactionEvolver < FieldImage > Reaction;
      Reaction reaction( epsilon );
      */

      // Explicit Reaction Evolver (with possible volume conservation)
      typedef ExplicitReactionEvolver<
          ImageContainerBySTLVector<Domain,double>,
          ImageContainerBySTLVector<Domain,double>
        > Reaction;
      ImageContainerBySTLVector<Domain,double> a( d );
      std::fill(a.begin(), a.end(), 0.0 );
      Reaction reaction( epsilon, a, balloon, flagWithCstVol );

      // Lie splitting
      LieSplittingEvolver< Diffusion, Reaction > phaseEvolver(diffusion, reaction);

      // Multi phase-field
      MultiPhaseField< LabelImage, FieldImage, LieSplittingEvolver<Diffusion,Reaction> > evolver(*labelImage, phaseEvolver);

      DGtal::trace.beginBlock( "Deformation (multi phase field)" );

      // Initial state export
      std::stringstream s;
      s << outputFiles << setfill('0') << std::setw(4) << 0;
#if   DIMENSION == 2
      drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
      writePartition( *labelImage, s.str(), outputFormat );
#endif

      // VTK export
        {
          VTKLightWriter<Domain> vtk(s.str(), labelImage->domain());
          for (size_t j = 0; j < evolver.getNumPhase(); ++j)
            {
              stringstream s_phase;
              s_phase << "phi" << setfill('0') << std::setw(2) << j;
              vtk << s_phase.str() << evolver.getPhase(j);
            }
        }

      // Time integration
      double sumt = 0;
      for (unsigned int i = 1; i <= max_step; ++i)
        {
          DGtal::trace.info() << "iteration # " << i << std::endl;

          // Update
          evolver.update( tstep );

          // Display
          if ( (i % disp_step) == 0 )
            {
              std::stringstream s;
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);
#if   DIMENSION == 2
	            drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
              writePartition( *labelImage, s.str(), outputFormat );
#endif

              // VTK export
              VTKLightWriter<Domain> vtk(s.str(), labelImage->domain());
              for (size_t j = 0; j < evolver.getNumPhase(); ++j)
                {
                  stringstream s_phase;
                  s_phase << "phi" << setfill('0') << std::setw(2) << j;
                  vtk << s_phase.str() << evolver.getPhase(j);
                }

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

          // Volume of each phase
          DGtal::trace.info() << ( dimension == 2 ? "Area: " : "Volume: " );
          for (size_t j = 0; j < evolver.getNumPhase(); ++j)
            {
              DGtal::trace.info() << "V(" << j << ") = " << getVolume<double>(evolver.getPhase(j));
            }
          DGtal::trace.info() << std::endl;

          DGtal::trace.info() << "Time spent: " << sumt << std::endl;
        }

      DGtal::trace.endBlock();

#if DIMENSION == 3
      // Interactive display after the evolution
      if ( flagWithVisu )
        displayPartition( *labelImage );
#endif

    }
  else if ( algo == "massiveMultiPhaseField" )
    {
      // Options's validity
      if (epsilon <= 0)
        {
          trace.error() << "epsilon should be greater than 0" << std::endl;
          return 1;
        }

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
      epsilon /= labelImage->extent()[0];
      tstep = epsilon*epsilon;

      MultiPhaseField2< LabelImage, FieldImage, ApproximatedMultiImage > evolver(*labelImage, epsilon, !no_dist);

      DGtal::trace.beginBlock( "Deformation (massive multi phase field)" );

      // Initial state export
      std::stringstream s;
      s << outputFiles << setfill('0') << std::setw(4) << 0;
#if   DIMENSION == 2
      drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
      writePartition( *labelImage, s.str(), outputFormat );
#endif

      // VTK export
        {
          VTKLightWriter<Domain> vtk(s.str(), labelImage->domain(),
              evolver.myRealExtent / ( labelImage->domain().upperBound() - labelImage->domain().lowerBound() + Point::diagonal() ) );
          /*
          for (size_t j = 0; j < evolver.getNumPhase(); ++j)
            {
              stringstream s_phase;
              s_phase << "phi" << setfill('0') << std::setw(2) << j;
              vtk << s_phase.str() << evolver.getPhase(j);
            }
          */
          vtk << "label" << *labelImage;
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
              const std::size_t label_cnt = evolver.updateLabels();

              // Display phase field informations
              evolver.dispInfos();

              // Export
              trace.beginBlock("Export");

              auto const implicit = makeFunctorConstImage( labelImage->domain(),
                [&evolver, epsilon] ( Point const& aPoint ) -> real
                  {
                    real max1 = 0, max2 = 0;
                    for ( auto value : evolver.getPhasesContainer()(aPoint) )
                      {
                        if ( value.second >= max1 ) { max2 = max1; max1 = value.second; }
                        else if ( value.second > max2 ) { max2 = value.second; }
                      }
                    //return max1 - max2;
                    return 2 * epsilon * std::atanh( std::min(max1 - max2, real(1)-1e-8) );
                  }
              );

              std::stringstream s;
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);
#if   DIMENSION == 2
              drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
              writePartition( *labelImage, s.str(), outputFormat );
#endif

              RawWriter< decltype(implicit) >::exportRaw<real>( s.str()+".imp.raw", implicit );
              RawWriter< LabelImage >::exportRaw<unsigned short int>( s.str()+".lab.raw", *labelImage );

              // VTK export
              VTKLightWriter<Domain> vtk(s.str(), labelImage->domain(),
                  evolver.myRealExtent / ( labelImage->domain().upperBound() - labelImage->domain().lowerBound() + Point::diagonal() ) );
              /*
              for (size_t j = 0; j < evolver.getNumPhase(); ++j)
                {
                  stringstream s_phase;
                  s_phase << "phi" << setfill('0') << std::setw(2) << j;
                  vtk << s_phase.str() << evolver.getPhase(j);
                }
              */
              vtk << "label" << *labelImage;
              //vtk << "implicit" << implicit;
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

              trace.endBlock();

              if ( label_cnt <= 0.000005 * disp_step * labelImage->domain().size() ) 
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
            << " -d " << labelImage->extent()[0]
            << " -d " << labelImage->extent()[1]
            << " -d " << labelImage->extent()[2]
            << " -S " << evolver.myRealExtent[0]
            << " -S " << evolver.myRealExtent[1]
            << " -S " << evolver.myRealExtent[2]
            //<< " -l " << ( s.str() + ".lab.raw" )
            //<< " -v no"
            << std::endl;
        }

#if DIMENSION == 3
      // Interactive display after the evolution
      if ( flagWithVisu )
        displayPartition( *labelImage );
#endif
    }
  else if ( algo == "localLevelSet" )
    {
      // Space
      KSpace ks;
      ks.init( d.lowerBound(), d.upperBound(), true );

      // Distance image...
      typedef ImageContainerBySTLVector<Domain,double> DistanceImage;
      // ... and extern data.
      DistanceImage g( d );
      std::fill(g.begin(), g.end(), 1.0);

      // Topological predicate
      typedef SimplePointHelper<LabelImage> TopologicalPredicate;
      TopologicalPredicate topologicalPredicate(*labelImage);

      // Frontier evolver
      DGtal::trace.beginBlock("Partition construction");
      PartitionEvolver<KSpace, LabelImage, DistanceImage, DistanceImage, TopologicalPredicate>
        e(ks, *labelImage, g, topologicalPredicate);

      DGtal::trace.info() << e << std::endl;
      DGtal::trace.endBlock();

      DGtal::trace.beginBlock( "Deformation (narrow band with topological control)" );

      // Initial state export
      std::stringstream s;
      s << outputFiles << setfill('0') << std::setw(4) << 0;
#if   DIMENSION == 2
      drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
      writePartition( *labelImage, s.str(), outputFormat );
#endif

      // Time integration
      double sumt = 0;
      for (unsigned int i = 1; i <= max_step; ++i)
        {
          DGtal::trace.info() << "iteration # " << i << std::endl;

          // Update
          e.update(tstep);

          // Display
          if ((i % disp_step) == 0)
            {
              std::stringstream s;
              s << outputFiles << setfill('0') << std::setw(4) << (i/disp_step);
#if   DIMENSION == 2
	            drawContours( *labelImage, s.str(), outputFormat );
#elif DIMENSION == 3
              writePartition( *labelImage, s.str(), outputFormat );
#endif
            }

          sumt += tstep;
          DGtal::trace.info() << "Time spent: " << sumt << std::endl;
        }

      DGtal::trace.endBlock();

#if DIMENSION == 3
      // Interactive display after the evolution
      if ( flagWithVisu ) displayPartition( *labelImage );
#endif

    }
  else
    {
      trace.error() << "unknown algo. Try option -h to see the available algorithms " << std::endl;
      return 1;
    }

  // Free
  delete( labelImage );

  return 0;
}

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
