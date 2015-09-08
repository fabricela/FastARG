/*
 * ARGnode.h
 *
 *  Created on: 2015-09-04
 *      Author: Éric Marcotte
 */

#ifndef ARGNODE_H_
#define ARGNODE_H_

#include <iostream>
#include <vector>

enum class ARGNodeType {LEAFSAMPLEMODE,MUTATIONODE,RECOMBINATIONODE,COALESCENSENODE};

class ARGnode {
public:

	static unsigned short int sizeOfSequence;
	static unsigned short int sizeOfSample;
	unsigned short int positionOfDataSource;
	unsigned short int secondDataPoint;
	ARGNodeType typeOfNode;
	std::vector<unsigned short int> mutationPositionList;
	std::vector<unsigned short int> ancestralMaterialPositionList;
	std::vector<std::vector<unsigned short int>> listOfLeafPositionForTIM; //2*SIZEOFSEQUENCE-2,std::vector<unsigned short int>(SIZEOFSAMPLE) );


	// Constructors
	// Constructor for a leaf node.
	ARGnode( std::vector<unsigned short int> mutationPositionList, unsigned short int positionOfTheNodeInData);

	// Constructor for a mutation node.
	ARGnode( unsigned short int positionOfSourceNode,
			 unsigned short int positionOfTheMutation , const std::vector<ARGnode>& ARGData );

	 // Constuctor for a recombination
	ARGnode( unsigned short int positionOfSourceNode , unsigned short int positionOfTheRecombination,
			 bool trueForleftSideIsAncestralMateriel , const std::vector<ARGnode>& ARGData );

	// Constructor for a coalescence node.
	ARGnode( unsigned short int positionOfSourceNode,
			const std::vector<ARGnode>& ARGData , unsigned short int positionOfTheSecondSouceNode);

	// Check for rule of 0,3,5.
	//Destructors
	virtual ~ARGnode()=default;

	//Functions
	std::vector<unsigned short int> possibleRecombinationPoint();

	//Overloading of << operator
	friend std::ostream& operator << (std::ostream &, const ARGnode&);
};

#endif /* ARGNODE_H_ */
