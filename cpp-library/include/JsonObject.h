#ifndef JsonObject_H
#define JsonObject_H

// todo: how to use 3rd party library in nested folders ?
//#include "rapidjson/document.h"
#include <rapidjson/document.h> // rapidjson
#include <rapidjson/filereadstream.h>

class JsonObject
{
public:
    JsonObject(std::string, int, bool);
    void setJson(std::string, int);
    // rapidjson::Document getJsonDoc(); // cannot return due to it being prohibited by document.h
    const char* getValue(std::string);

    bool getAssertionStatus();
    void setAssertionStatus(bool);
private:
    rapidjson::Document jsonDoc;
    bool checkAssertions;
};

#endif