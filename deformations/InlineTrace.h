#pragma once

#include <ostream>
#include <string>
#include <sstream>
#include <DGtal/base/Clock.h>

namespace DGtal
{

class InlineTrace
{
public:
  InlineTrace( std::ostream& aStream )
      : myStream(aStream)
    {}

  void beginBlock( std::string const& aDesc = "Block" )
    {
      myDesc = aDesc;
      mySteps.clear();
      dispMe();
      myGlobalClock.startClock();

    }

  double endBlock()
    {
      myTime = myGlobalClock.stopClock();

      if ( isRunningStep() )
        endStep();

      double steps_time = 0;
      for ( auto const& step : mySteps )
        steps_time += step.second;
      mySteps.emplace_back( Step{ "Other", myTime - steps_time } );

      dispMe(true);
    }

  void beginStep( std::string const& aShortDesc )
    {
      if ( isRunningStep() ) 
        endStep();

      if ( mySteps.size() > 0 ) myStream << " ; ";
      myStream << aShortDesc << "...\r" << std::flush;

      mySteps.emplace_back( Step{aShortDesc, -1} );
      myStepClock.startClock();
    }

  double endStep()
    {
      mySteps.back().second = myStepClock.stopClock();
      dispMe(false);
    }

private:
  using Step = std::pair<std::string, double>;

  inline
  bool isRunningStep()
    {
      return mySteps.size() > 0 && mySteps.back().second == -1;
    }

  std::string formatTime( double aTime )
    {
      std::stringstream format;
      format << setprecision(0) << std::fixed << aTime << "ms";
      return format.str();
    }

  void dispStep( Step const& aStep )
    {
      myStream << aStep.first << ":" << formatTime(aStep.second);
    }

  void dispMe( bool isLast = false )
    {
      if ( isLast )
          myStream << "\r[" << myDesc << ":" << formatTime(myTime) << "] ";
      else
          myStream << "[" << myDesc << "] ";

      bool is_first = true;
      for ( auto const& step : mySteps )
        {
          myStream << (is_first ? "" : " ; ") << step.first << ":" << formatTime(step.second);
          is_first = false;
        }

      if ( isLast )
        myStream << std::endl;
      else
        myStream << std::flush;
    }

  std::ostream& myStream;
  std::string myDesc;
  std::list<Step> mySteps;
  Clock myGlobalClock, myStepClock;
  double myTime;

};

  InlineTrace itrace(std::cout);

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
