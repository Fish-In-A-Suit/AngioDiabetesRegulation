// SequenceComparator.h

#ifndef SequenceComparator_H
#define SequenceComparator_H

#include <vector>
#include <string>
#include <iostream>

class SequenceComparator {
	std::vector<std::string> miRNAsequences;
	std::vector<std::string> mRNAsequences;

public:
	SequenceComparator(std::string, std::string);
	~SequenceComparator();
	void load_miRNA_sequences(std::string);
};

#endif