// Logger.cpp

#include "Logger.h"

/**
 * All inline static functions should be declared AND defined in .h file. Defining in .cpp file raises "function used but never defined" error.
 * 
void Logger::setLevel(Constants::LogLevels level) {
    logLevel = level;
}

Constants::LogLevels Logger::getLevel() {
    return logLevel;
}

void Logger::setName(std::string name) {
    loggerName = name;
}

std::string Logger::getName() {
    return loggerName;
}

void Logger::info(std::string message) {
    if ((logLevel == Constants::LogLevels::INFO) || (logLevel == Constants::LogLevels::DEBUG)) {
        std::cout << message << std::endl;
    }
}

void Logger::debug(std::string message) {
    if (logLevel == Constants::DEBUG) {
        std::cout << message << std::endl;
    }
}

void checkType(auto *variable) {
    checkType(variable, "");
}

void checkType(auto *variable, std::string marker) {
    std::stringstream sstream;
    sstream << marker << " type: " << typeid(variable).name() << "";
    Logger::debug(sstream.str());
}
*/