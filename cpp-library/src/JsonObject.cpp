#include "JsonObject.h"


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

bool JsonObject::getAssertionStatus() {
    return checkAssertions;
}

void JsonObject::setAssertionStatus(bool newStatus) {
    checkAssertions = newStatus;
}