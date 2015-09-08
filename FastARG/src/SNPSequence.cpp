/*
 * SNPSequence.cpp
 *
 *  Created on: 2015-05-20
 *      Author: Éric Marcotte
 */


#include "SNPSequence.h"

//Constructors
//SNPSequence::SNPSequence() {}

SNPSequence::SNPSequence(int ID , bool asTIM, bool isACase) {
	this->sampleId = ID;
	this->asTIM = asTIM;
	this->isACase = isACase;
}

SNPSequence::SNPSequence(int ID , bool asTIM, bool isACase,std::vector<unsigned short int> mutationsPositions) {
	this->mutationPositionList = mutationsPositions;
	this->sampleId = ID;
	this->asTIM = asTIM;
	this->isACase = isACase;
}


//Destructors
SNPSequence::~SNPSequence() {}

//Operators overlaoding
std::ostream& operator << (std::ostream & out, const SNPSequence & aSeq){
	out << "Id:" << aSeq.sampleId
			<< " Case: " << aSeq.asTIM
			<< " control:" << aSeq.isACase
			<< " Total number of mutation: "
			<< aSeq.mutationPositionList.size() << " mutation position: ";

	for(std::vector<unsigned short int>::const_iterator it = aSeq.mutationPositionList.begin();
			it != aSeq.mutationPositionList.end();
			++it){
		out << *it << ' ';
	}
	return out;
}

bool SNPSequence::operator == (const SNPSequence & ASeq) const{
	return this->mutationPositionList == ASeq.mutationPositionList;
}

bool operator <(const SNPSequence& aSeq, const SNPSequence& aSeq2 ){
	return  aSeq.mutationPositionList.size() < aSeq2.mutationPositionList.size();
}
