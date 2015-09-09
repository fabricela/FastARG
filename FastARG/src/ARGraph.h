/*
 * ARGraph.h
 *
 *  Created on: 2015-09-05
 *      Author: Éric Marcotte
 */

#ifndef ARGRAPH_H_
#define ARGRAPH_H_
#include "ARGnode.h"
#include "SNPSequence.h"

class ARGraph {
public:
	unsigned short int sizeOfSequence;
	unsigned short int sizeOfSample;

	std::vector<SNPSequence> dataSequences;
	std::vector<ARGnode> ARGnodes;
	unsigned short int distanceBetweenMarquers;

	std::vector<unsigned short int> positionOfLiveNodes;
	std::vector<unsigned short int> leafPositionWithTheTIM;
	std::vector<unsigned short int> sampleIDThatAreCases;
	std::vector<int> statistics;


	// Constructors
	ARGraph(const std::vector<SNPSequence>& dataSequences,unsigned short int,unsigned short int);

	// Destructors
	// Check for rule of 0,3,5.
	virtual ~ARGraph() = default;

	// ARG generaters
	void constructExampleARG();
	void randomArgGenerator();

	// Functions
	void generateACoalescentEvent(unsigned short int, unsigned short int);
	void generateARecombinationEvent(unsigned short int, unsigned short int);
	void generetaAMutationEvent(unsigned short int, unsigned short int);
	void generateStatisticData(); // Specialized for dominant or recessive gene.
	void nodeTests();
	bool generateARandomCoalescentEvent();
	bool generateARandomMutationEvent();
	bool generateARandomRecombinationEvent();
	void resetARG();
	void printLivesNodes();
	void printARGSummary();
	void printEvents();
};

#endif /* ARGRAPH_H_ */
