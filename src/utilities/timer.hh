#ifndef DUNE_DYCAP_TIMER_HH
#define DUNE_DYCAP_TIMER_HH

#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

namespace Dune {

  namespace Dycap{

    /*! \file
      \brief A simple timing class.
    */

      /** \brief A simple stop watch

        This class reports the elapsed user-time, i.e. time spent computing,
        after the last call to Timer::reset(). The results are seconds and
        fractional seconds. Note that the resolution of the timing depends
        on your OS kernel which should be somewhere in the milisecond range.

        The class is basically a wrapper for the libc-function getrusage()

        Taken from the DUNE project www.dune-project.org

    */
    class DycapTimer
    {
    public:
      //! A new timer, start immediately
      DycapTimer () throw(TimerError)
      {
        reset();
      }

      //! Reset timer
      void reset() throw (TimerError)
      {
        cstart = std::chrono::high_resolution_clock::now();
      }

      //! Get elapsed user-time in seconds
      double elapsed () const throw (TimerError)
      {
	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - cstart);
	return time_span.count();
      }

    private:
      std::chrono::high_resolution_clock::time_point cstart;
    }; // end class Timer

  } // end namespace
}

#endif
