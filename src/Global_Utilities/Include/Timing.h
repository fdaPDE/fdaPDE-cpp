#ifndef __TIMING_H__
#define __TIMING_H__

#include <iostream>
#include <iomanip>
#include <ctime>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

//! a class to measure the code performance in terms of time
class timer {
public:
  void start() {
    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    begin.tv_sec = mts.tv_sec;
    begin.tv_nsec = mts.tv_nsec;

    #else
    clock_gettime(CLOCK_REALTIME, &begin);
    #endif
  }

  timespec stop() {
    timespec end;

    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    end.tv_sec = mts.tv_sec;
    end.tv_nsec = mts.tv_nsec;

    #else
    clock_gettime(CLOCK_REALTIME, &end);
    #endif

    timespec difference = diff(begin, end);
    //Rprintf("It took %u.%09us\n", difference.tv_sec, difference.tv_nsec);
//    std::cout << "It took " << difference.tv_sec << "."
//              << std::setfill('0') << std::setw(9) << difference.tv_nsec
//              << "s" << std::endl;
    return difference;
  }
private:
  timespec diff(timespec start, timespec end)
  {
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
      temp.tv_sec = end.tv_sec-start.tv_sec-1;
      temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
      temp.tv_sec = end.tv_sec-start.tv_sec;
      temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
  }

private:
  timespec begin;
};

#endif
