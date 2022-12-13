#include "test-cpp-lib-hello.h"
#include "StringUtils.h"
#include "HighResolutionTimeManager.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <chrono>
#include <json/json.h> //jsoncpp
#include <rapidjson/document.h> //rapidjson
#include "rapidjson/filereadstream.h"
#include <cstdio>

// #include <test-cpp-library/include/test-cpp-lib-hello.h>

using namespace std;

/*
Workflow for miRNA prediction:
  (1) input parameters:
        - productNames: the uniprot ids of products
        - mRNAs: array of mRNA molecules to be analysed (each corresponds to one productName)
        - productScores: the scores of the products
        - miRNALength: the length of the miRNAs that will be computed. Immensely increases required computation power after 14+
        - thresholdToAccept:
        - targetFolder:
        - debugLoopBreak:

  (2) Given miRNALength, generate all possible (unique) mRNA substrings (from all mRNAs) and store them in mRNA_unique_substrings
      Hash the substrings into bits --> mRNA_unique_substrings_hashed. If this is not done, you would have to do complex bit-stepping in part (3)
      (you'd have to step through the mRNA sequence bit by bit and append miRNA at different starting bits, first in a forward direction and then
      backwards too)

  miRNA_mRNA_fitting_dict =
  {
    {
        "miRNA": "xxxxx",
        "mRNAsubstrings": ["xxxxxx", "aaaaaa", "bbbbbb"]
        "mRNAsubstringsScores": [1.0, 0.25, 0.32]
    }
    {
        "miRNA": "xxxxx",
        "mRNAsubstrings": ["xxxxxx", "aaaaaa", "bbbbbb"]
        "mRNAsubstringsScores": [1.0, 0.25, 0.32]
    }
  }

  numberOfPermutations = compute total number of permutations (mathematically)
  int i = 0
  // also compute timing (averageNestedForLoopTime = totalNestedForLoopTime/i) and then find ETA from AverageNestedForLoopTime * (numberOfPermutations - i)
  (3) top level for loop: loop through possible miRNA permutations of ["A", "T", "C", "G"] of given miRNALength
      for each miRNA permutation (in a bit-format):
        matchedSubstrings = []
        matchedScores = [] // the scores of the matchedSubstrings
        for each mRNA substring:
            score = check number of matching nucleotides (if reading from left to right, 2 bits have to be the same consecutively for a single nucleotide match)
            if score > threshold:
                matchedSubstrings.append(mRNA)
                matchedScores.append(score)
        append miRNA, mRNAsubstrings and mRNAsubstringScores into miRNA_mRNA_fitting_dict
        i++
      order & save miRNA_mRNA_fitting_dict
*/

// todo: port into OOP style
void predict_miRNA_values(list<string> t)
{
}

long test_sequence_container_speed()
{
    // Use auto keyword to avoid typing long type definitions
    auto start = std::chrono::high_resolution_clock::now();
    int i = 111111111;

    // todo: test vectors, lists, forward lists and deques for speed comparisons
}

int main()
{
    // "functional-programming" way of computing elapsed time:
    // auto startTime = high_resolution_clock::now();
    // cout << compute_elapsed_time(startTime, Constants::TimeUnits::NANOSECONDS) << endl;
    std::string helloJim = generateHelloString("Aljosa & Ladi");
    std::cout << helloJim << std::endl;

    vector<string> msg{"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};
    int i = 0;
    for (const string &word : msg)
    {
        cout << word << " ";
        ++i;
    }
    cout << endl;

    // *** JsonCpp json parsing ***:
    // this doesnt work with a relative filepath, but absolute filepath is ok.
    // BUG: JSON::READER DOESN'T HANDLE RELATIVE FILEPATHS!
    //std::ifstream json_file_ifs("C:\\Aljosa\\Development\\Unity-Github\\AngioDiabetesRegulation\\src_data_files\\test.json", std::ifstream::binary);
    //Json::Reader reader;
    //Json::Value root;
    //reader.parse(json_file_ifs, root, false); // // reader can also read strings
    //cout << root << endl;
    //cout << "Book: " << root["book"].asString() << endl;
    //std::string encoding = root.get("encoding", "UTF-8").asString();
    //std::cout << encoding << "\n";

    // *** RapidJSON parsing ***
    HighResolutionTimeManager hrtm;
    FILE *fp = fopen("C:\\Aljosa\\Development\\Unity-Github\\AngioDiabetesRegulation\\src_data_files\\test.json", "rb");
    char readBuffer[65536];
    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    rapidjson::Document document;
    document.ParseStream(is);
    fclose(fp);
    cout << "Took " << hrtm.getElapsedTime(Constants::MICROSECONDS) << " to parse file." << endl;

    assert(document.HasMember("book"));
    assert(document["book"].IsString());
    printf("book = %s\n", document["book"].GetString());
}

/**
 * A "functional-programming" way of coding to compute elapsed time -> ported to OOP with the HighResolutionTimeManager class.
 *
 * Computes elapsed program runtime by subtracting current time (computed at this function call) from the supplied start time,
 * returning the amount of units (seconds, microseconds, etc.) passed.
 *
 * @param start: the start time computed by std::chrono::high_resolution_clock::now()
 * @param unit: either "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes, "hours" from Constants.TimeUnits
 *
long compute_elapsed_time(std::chrono::high_resolution_clock::time_point start, Constants::TimeUnits timeUnit)
{
    auto stop = std::chrono::high_resolution_clock::now();
    int i = 4;
    auto duration = std::chrono::duration_cast<nanoseconds>(stop - start); // initialise auto var
    // namespaces used for code portability
    // Constants::TimeUnits timeUnit = timeUnit;
    switch (timeUnit)
    {
    case Constants::TimeUnits::NANOSECONDS:
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        break;
    case Constants::TimeUnits::MICROSECONDS:
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        break;
    case Constants::TimeUnits::MILLISECONDS:
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        break;
    case Constants::TimeUnits::SECONDS:
        duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        break;
    case Constants::TimeUnits::HOURS:
        duration = std::chrono::duration_cast<std::chrono::hours>(stop - start);
        break;
    }
    return duration.count();
}
*/