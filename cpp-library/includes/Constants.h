// Constants.h

#ifndef Constants_H
#define Constants_H

class Constants {
public:
    // reasons for implementing enum in .h: https://stackoverflow.com/questions/1284529/enums-can-they-do-in-h-or-must-stay-in-cpp
    enum TimeUnits
    {
        NANOSECONDS,
        MICROSECONDS,
        MILLISECONDS,
        SECONDS,
        MINUTES,
        HOURS
    };

    enum LogLevels {
        DEBUG,
        INFO
    };
};

#endif