// HighResolutionTimeManagerV2.cpp

#include "HighResolutionTimeManagerV2.h"

HighResolutionTimeManagerV2::HighResolutionTimeManagerV2() {
    setStartTime();
}

void HighResolutionTimeManagerV2::setStartTime() {
    startTime = std::chrono::high_resolution_clock::now();
}

std::string HighResolutionTimeManagerV2::getElapsedTime(Constants::TimeUnits timeUnit, bool print) {
    auto stopTime = std::chrono::high_resolution_clock::now();
    double duration;
    double duration_sec;
    std::ostringstream resultStrStream;

    switch(timeUnit) {
        case Constants::TimeUnits::NANOSECONDS:
            duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stopTime - startTime).count();
            duration_sec = duration * 1e-9;
            resultStrStream << "Time: " << std::fixed << duration_sec << std::setprecision(9) << std::setprecision(1) << "s (" << duration << " ns)";
            break;
        case Constants::TimeUnits::MICROSECONDS:
            duration = std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime).count();
            duration_sec = duration * 1e-6;
            resultStrStream << "Time: " << std::fixed << duration_sec << std::setprecision(6) << std::setprecision(1) << "s (" << duration << " mcrs)";
            break;
        case Constants::TimeUnits::MILLISECONDS:
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count();
            duration_sec = duration * 1e-3;
            resultStrStream << "Time: " << std::fixed << duration_sec << std::setprecision(6) << std::setprecision(2) << "s (" << duration << " ms)";
            break;
        case Constants::TimeUnits::SECONDS:
            duration = std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime).count();
            duration_sec = duration;
            resultStrStream << "Time: " << std::fixed << duration_sec << std::setprecision(6) << "s";
            break;
        case Constants::TimeUnits::HOURS:
            duration = std::chrono::duration_cast<std::chrono::hours>(stopTime - startTime).count();
            duration_sec = duration * 3600;
            resultStrStream << "Time: " << std::fixed << duration << std::setprecision(9) << std::setprecision(1) << "hours (" << duration_sec << " s)";
            break;
    }
    if (print) {
        std::cout << resultStrStream.str() << std::endl;
    } else {
        return resultStrStream.str();
    }
}