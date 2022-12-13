// HighResolutionTimeManager.h

#ifndef HighResolutionTimeManager_H
#define HighResolutionTimeManager_H

#include <chrono>
#include <cmath>
#include "StringUtils.h"
#include "Constants.h"

// never use 'using namespace <...>' in .h !!!

// note: DO NOT USE () next to class defition here !
class HighResolutionTimeManager {
    std::chrono::high_resolution_clock::time_point startTime;
    StringUtils StringUtils;

public:
    HighResolutionTimeManager();
    void setStartTime();
    long getElapsedTime(Constants::TimeUnits, bool);
    std::string getElapsedTime(Constants::TimeUnits);

private:
    std::string formatTimeValue(long, Constants::TimeUnits);
    std::string _divideTimeAndReturn(long, long, int, std::string, std::string, int);
};

#endif