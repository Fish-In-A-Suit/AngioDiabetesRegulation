// Logger.h

#ifndef Logger_H
#define Logger_H

#include "Constants.h"
#include <string>
#include <iostream>
#include <typeinfo>
#include <sstream>
#include <rapidjson/document.h>

class Logger {
    inline static Constants::LogLevels logLevel;
    inline static std::string loggerName;

public:
    inline static void setLevel(Constants::LogLevels level) {
        logLevel = level;
    };

    inline static Constants::LogLevels getLevel() {
        return logLevel;
    };

    inline static void setName(std::string name) {
        loggerName = name;
    };

    inline static std::string getName() {
        return loggerName;
    };

    inline static void info(std::string message) {
        if ((logLevel == Constants::LogLevels::INFO) || (logLevel == Constants::LogLevels::DEBUG)) {
            std::cout << message << std::endl;
        }
    };

    inline static void debug(std::string message) {
       if (logLevel == Constants::DEBUG) {
            std::cout << message << std::endl;
        } 
    };

    inline static void checkType(auto *variable) {
        checkType(variable, "");
    };

    inline static void checkType(auto *variable, std::string marker) {
        std::stringstream sstream;
        if (marker == "") {
            sstream << "Type: " << typeid(variable).name() << "";
        } else {
            sstream << marker << " type: " << typeid(variable).name() << "";
        }
        
        Logger::debug(sstream.str());
    };

private:
    Logger();
};

#endif