//
//  debug.h
//  cesna
//
//  Created by 金宇超 on 14-12-27.
//  Copyright (c) 2014年 ccjinyc. All rights reserved.
//
#ifndef __cesna__debug__
#define __cesna__debug__
#include <iostream>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;


#define DEBUG 1
#define USE_PROFILING 1
#define USE_TIMESTAMP 1
#if USE_PROFILING
#define PROFILE(msg,code)	{\
using namespace std::chrono;\
auto t1 = Clock::now();\
{ code }\
auto t2 = Clock::now();\
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);\
std::cout << msg << " ";\
std::cout << "DONE " << time_span.count() << " seconds." << std::endl;\
}
#else
#define PROFILE(msg,code) code
#endif //USE_PROFILING


#if USE_TIMESTAMP
#define SETUPSTAMP \
using namespace std::chrono;\
auto t1 = Clock::now();
#define STAMP(msg, code) {\
using namespace std::chrono;\
{code}\
auto t2 = Clock::now();\
duration<double> time_span = duration_cast<duration<double>>(t2 - t1);\
std::cout << msg << " " << time_span.count() << " seconds." << std::endl;\
}
#else
#define STAMP(msg, code) code
#define SETUPSTAMP {}
#endif //USE_TIMESTAMP

#endif // __cesna__debug__
