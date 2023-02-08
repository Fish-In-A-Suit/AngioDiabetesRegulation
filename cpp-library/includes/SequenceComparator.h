// SequenceComparator.h

#ifndef SequenceComparator_H
#define SequenceComparator_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "StringUtils.h"

class SequenceComparator {
	std::vector<std::string> miRNA_sequences; // contains all of the miRNA sequences
	std::vector<std::string> miRNA_miRDB_ids; // contains all of the respective miRNA miRDB tags (eg MI000001)
	std::vector<std::string> miRNA_names; // contains all of the respective miRNA ids (eg. hsa-let-7a-1)
	std::vector<std::string> mRNA_sequences;

public:
	SequenceComparator(std::string, std::string);
	~SequenceComparator();
	void load_miRNA_sequences(std::string);

private:
	std::vector<std::vector<std::string>> process_miRNAsequences_file(std::string); // processes mirbase_miRNA_hsa-only.txt
	std::vector<std::string> process_mRNAsequences_file(std::string);
};

#endif