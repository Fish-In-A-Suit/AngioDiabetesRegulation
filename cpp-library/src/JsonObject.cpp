#include "JsonObject.h"


/**
 * Creates a json object.
 * @param filepath: the filepath to the json file
 * @param readBufferSize: the size of the read buffer
 * @param checkAsserts: if True, assertions will be checked when querying values (slower but safer)
*/
JsonObject::JsonObject(std::string filepath, int readBufferSize, bool checkAsserts)
{
    setJson(filepath, readBufferSize);
    checkAssertions = checkAsserts;
}

void JsonObject::setJson(std::string filepath, int readBufferSize) {
    FILE *fp = fopen(filepath.c_str(), "rb");
    char readBuffer[readBufferSize];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    jsonDoc.ParseStream(is);
    fclose(fp);
}

const char* JsonObject::getValue(std::string key){
    const char* keyCstr = key.c_str();

    if (checkAssertions) {
        // assertions = true, do checking
        assert(jsonDoc.HasMember(keyCstr));
        assert(jsonDoc[keyCstr].IsString()); // todo: maybe remove this
        return jsonDoc[keyCstr].GetString();
    } else {
        // assertions = false, do no checking
        return jsonDoc[keyCstr].GetString();
    }
}

/**
 * Converts private variable jsonDoc to string format.
 * 
 * @param keepIndentation: If true, indents are preserved in the resulting string. If false, returns a compat json structure without indentation.
 * !!! Note that enabling indentation causes a significant performance drop with bigger documents (0.196s without indents -> 2.25s with indents) !!!
*/
std::string JsonObject::toString(bool keepIndentation) {
    rapidjson::StringBuffer buffer;
    std::string json_string;

    if (keepIndentation) {
        // keep indents with pretty writer
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(buffer);
        jsonDoc.Accept(writer);
    } else {
        // compact with regular writer
        rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
        jsonDoc.Accept(writer);
    }
    json_string = buffer.GetString();
    return json_string;
}

bool JsonObject::getAssertionStatus() {
    return checkAssertions;
}

void JsonObject::setAssertionStatus(bool newStatus) {
    checkAssertions = newStatus;
}