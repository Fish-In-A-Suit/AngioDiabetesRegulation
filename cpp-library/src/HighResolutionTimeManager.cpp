// HighResolutionTimeManager.cpp

#include "HighResolutionTimeManager.h"
#include <iostream>

HighResolutionTimeManager::HighResolutionTimeManager() {
    setStartTime();
}

void HighResolutionTimeManager::setStartTime() {
    startTime = std::chrono::high_resolution_clock::now();
}

/**
 * Computes elapsed program runtime by subtracting current time (computed at this function call) from the supplied start time,
 * returning the amount of units (seconds, microseconds, etc.) passed.
 *
 * @param timeUnit: either "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes, "hours" from Constants.TimeUnits
 * @param _o_rawLong: a function overloading parameter: servers no algorithmical value except overloading (and calling)
 * long getElapsedTime from the string getElapsedTime() function
 * 
 * @returns: Elapsed time in specified time units as long.
 */
long HighResolutionTimeManager::getElapsedTime(Constants::TimeUnits timeUnit, bool _o_rawLong) {
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - startTime); // initialise auto var
    
    switch (timeUnit)
    {
    case Constants::TimeUnits::NANOSECONDS:
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - startTime);
        break;
    case Constants::TimeUnits::MICROSECONDS:
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - startTime);
        break;
    case Constants::TimeUnits::MILLISECONDS:
        // cout << "td duration = ", duration.count() << endl;
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - startTime);
        break;
    case Constants::TimeUnits::SECONDS:
        duration = std::chrono::duration_cast<std::chrono::seconds>(stop - startTime);
        break;
    case Constants::TimeUnits::HOURS:
        duration = std::chrono::duration_cast<std::chrono::hours>(stop - startTime);
        break;
    }

    long long_duration = (long)duration.count();
    std::cout << duration.count() << std::endl;
    return long_duration;
}

/**
 * Computes elapsed program runtime by subtracting current time (computed at this function call) from the supplied start time,
 * returning the amount of units (seconds, microseconds, etc.) passed.
 *
 * @param timeUnit: either "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes, "hours" from Constants.TimeUnits
 *
 * @returns: Elapsed time in specified time in a user-friendly string notation.
 */
std::string HighResolutionTimeManager::getElapsedTime(Constants::TimeUnits timeUnit) {
    return formatTimeValue(getElapsedTime(timeUnit, true), timeUnit);
}

std::string HighResolutionTimeManager::formatTimeValue(long duration, Constants::TimeUnits timeUnit) {
    long divider = 1000L;
    float divided;
    int upflowThreshold = 1;
    int decimalAccuracy = 2;
    
    if (timeUnit == Constants::SECONDS)
    {
        divider = 60;
    }
    if (timeUnit == Constants::HOURS)
    {
        divider = 24;
    }

    switch (timeUnit)
    {
    case Constants::TimeUnits::NANOSECONDS:
        // if more than 1000, format to mcrs
        return _divideTimeAndReturn(duration, divider, upflowThreshold, "mcrs", "ns", decimalAccuracy);
        break;
    case Constants::TimeUnits::MICROSECONDS:
        // if more than 1000, format to ms (0.01xxx ms, 0.1xxx ms, 1.xxx ms)
        return _divideTimeAndReturn(duration, divider, upflowThreshold, "ms", "mcrs", decimalAccuracy);
        break;
    case Constants::TimeUnits::MILLISECONDS:
        // if more than 1000, format to seconds
        return _divideTimeAndReturn(duration, divider, upflowThreshold, "s", "ms", decimalAccuracy);
        break;
    case Constants::TimeUnits::SECONDS:
        // if more than 60, format to minutes, if more than 3600 format to hours
        return _divideTimeAndReturn(duration, divider, upflowThreshold, "h", "s", decimalAccuracy);
        break;
    case Constants::TimeUnits::HOURS:
        // if more than 24, format to days
        return _divideTimeAndReturn(duration, divider, upflowThreshold, "d", "h", decimalAccuracy);
        break;
    }
    return "-1";
}

/**
 * Example: You have 14.500ns (duration), use 1000 as divider. First, divide duration with divider (to get divided), in our
 * example divided = 14,5. If divided is greater than upflowThreshold (in our example: 1), then it will "upflow" (or "overflow") into
 * the next number ie 14,5 micro seconds (upflowMarker "mcrs" used). If divided is lower than upflowThreshold, then the final result
 * stays in the same numerical category and downflowMarker (in our example: "ns") is used.
 *
 * Using 1000 as divider and 1 as upflowThreshold means: all values lower than 1000 will stay in the same numerical category (ie 976ns),
 * whereas all values higher than 1000 will move into a higher numerical category (ie. 1050ns -> 1,05mcrs)
 * 
 * decimalAccuracy specifies how many decimals are used in the return result
 */
std::string HighResolutionTimeManager::_divideTimeAndReturn(long duration, long divider, int upflowThreshold, std::string upflowMarker, std::string downflowMarker, int decimalAccuracy) {
    double divided = duration / divider;
    divided = std::ceil(divided * 100) / 100;
    if (divided > upflowThreshold)
    {
        // problem: std::to_string doesn't allow to set precision (always xx.000000 decimals, I want only 2 decs precision)
        // return std::to_string(divided) + " " + upflowMarker;
        
        // changed conversion, to_string is now static
        //return StringUtils.to_string(divided, decimalAccuracy) + " " + upflowMarker;
        // if it upflows, after first division, another division is necessarry
        return StringUtils::to_string(divided/divider, decimalAccuracy) + " " + upflowMarker;
    }
    else
    {
        // return std::to_string(duration) + " " + downflowMarker;
        //return StringUtils.to_string(divided, decimalAccuracy) + " " + downflowMarker;
        return StringUtils::to_string(divided, decimalAccuracy) + " " + downflowMarker;
    }
}