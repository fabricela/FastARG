/*
 * SNPSequence.h
 *
 *  Created on: 2015-05-20
 *      Author: Éric Marcotte
 */

#ifndef SNPSEQUENCE_H_
#define SNPSEQUENCE_H_

#include <iostream>

#include <vector>         // std::bitset

class SNPSequence {
public:
	std::vector<unsigned short int> mutationPositionList;
	int sampleId;
	bool asTIM;
	bool isACase;


	//Constructors
	//SNPSequence();
	SNPSequence(int,bool,bool);
	SNPSequence(int,bool,bool,std::vector<unsigned short int>);

	//Destructors
	virtual ~SNPSequence();

	//Overloading of << operator
	friend std::ostream& operator << (std::ostream &, const SNPSequence &);
	//Overloading of < operator
	friend bool operator < (const SNPSequence &,const SNPSequence &);
	//Overloading of == operator
    bool operator == (const SNPSequence &) const;

};

#endif /* SNPSEQUENCE_H_ */
