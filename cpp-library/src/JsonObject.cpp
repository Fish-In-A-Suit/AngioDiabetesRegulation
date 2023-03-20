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

JsonObject::~JsonObject() {
    // destructor implementation
    std::cout << "JsonObject for " << filepath << " called." << std::endl;
}

void JsonObject::setJson(std::string filepath, int readBufferSize) {
    this->filepath = filepath; 
    FILE *fp = fopen(filepath.c_str(), "rb");
    char* readBuffer = new char[readBufferSize];
    //char readBuffer[readBufferSize]; // this is flagged as an error, since stack-based array have to have their size known at compile time.
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    jsonDoc.ParseStream(is);
    fclose(fp);
    delete[] readBuffer; // delete the readBuffer after usage
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