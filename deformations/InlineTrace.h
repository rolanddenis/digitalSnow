#pragma once

#include <ostream>
#include <string>
#include <sstream>
#include <boost/timer/timer.hpp>

namespace DGtal
{

class InlineTrace
{
public:
  InlineTrace( std::ostream& aStream )
      : myStream(aStream)
    {
      myGlobalClock.stop();
      myStepClock.stop();
    }

  void beginBlock( std::string const& aDesc = "Block" )
    {
      myDesc = aDesc;
      mySteps.clear();
      dispMe();
      myGlobalClock.start();

    }

  double endBlock()
    {
      myGlobalClock.stop();
      myTime = myGlobalClock.elapsed();

      if ( isRunningStep() )
        endStep();

      Time other_time = myTime;
      for ( auto const& step : mySteps )
        {
          other_time.wall -= step.second.wall;
          other_time.user -= step.second.user;
          other_time.system -= step.second.system;
        }
      mySteps.emplace_back( Step{ "Other", other_time } );

      dispMe(true);

      return static_cast<double>(myTime.wall)/1e6;
    }

  void beginStep( std::string const& aShortDesc )
    {
      if ( isRunningStep() ) 
        endStep();

      if ( mySteps.size() > 0 ) myStream << " ; ";
      myStream << aShortDesc << "...\r" << std::flush;

      mySteps.emplace_back( Step{aShortDesc, {-1,-1,-1} } );
      myStepClock.start();
    }

  double endStep()
    {
      myStepClock.stop();
      mySteps.back().second = myStepClock.elapsed();
      dispMe(false);
      return static_cast<double>(mySteps.back().second.wall)/1e6;
    }

private:
  using Time = boost::timer::cpu_times;
  using Step = std::pair<std::string, Time>;

  inline
  bool isRunningStep()
    {
      return mySteps.size() > 0 && mySteps.back().second.wall == -1;
    }

  std::string formatTime( Time aTime )
    {
      std::stringstream format;
      format  << setprecision(0) << std::fixed << static_cast<double>(aTime.wall)/1e6 << "ms"
              << "(x" << setprecision(1) << std::fixed << static_cast<double>(aTime.user+aTime.system)/aTime.wall << ")";
      return format.str();
    }

  void dispStep( Step const& aStep )
    {
      myStream << aStep.first << ":" << formatTime(aStep.second);
    }

  void dispMe( bool isLast = false )
    {
      if ( isLast )
          myStream << "\r[\e[1m" << myDesc << "\e[0m:" << formatTime(myTime) << "] ";
      else
          myStream << "[\e[1m" << myDesc << "\e[0m] ";

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
  boost::timer::cpu_timer myGlobalClock, myStepClock;
  Time myTime;

};

  InlineTrace itrace(std::cout);

} // namespace DGtal

/* GNU coding style */
/* vim: set ts=2 sw=2 expandtab cindent cinoptions=>4,n-2,{2,^-2,:2,=2,g0,h2,p5,t0,+2,(0,u0,w1,m1 : */
