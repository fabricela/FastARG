/*
 * ARGnode.cpp
 *
 *  Created on: 2015-09-04
 *      Author: Éric Marcotte
 */

#define VALIDATION

#include "ARGnode.h"
#include <iostream>
#include <limits>
#include <algorithm>

// Constructors
/*********************************************************************************************/
/**                                   LEAF NODE CONSTRUCTOR                                 **/
/*********************************************************************************************/
ARGnode::ARGnode( std::vector<unsigned short int> mutationPositionList , unsigned short int positionOfTheNodeInData ) {

	std::vector<unsigned short int> vectorFor1PositionOfTheTim;

#ifdef VALIDATION
	if ( this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() ){
		std::cerr << "this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() in ARGnode::ARGnode";
		throw "this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() in ARGnode::ARGnode";
	}
	if (this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() ){
		std::cerr << "this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() in ARGnode::ARGnode";
		throw "this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() in ARGnode::ARGnode";
	}
#endif
	this->typeOfNode =  ARGNodeType::LEAFSAMPLEMODE;
	this->mutationPositionList = mutationPositionList;
	this->positionOfDataSource = positionOfTheNodeInData;
	this->secondDataPoint = std::numeric_limits<unsigned short int>::max();
	vectorFor1PositionOfTheTim.push_back(positionOfTheNodeInData);
	for (unsigned short int i = 0 ; i < 2*this->sizeOfSequence-2 ; ++i){
		this->listOfLeafPositionForTIM.push_back(vectorFor1PositionOfTheTim);
	}
}

/*********************************************************************************************/
/**                               MUTATION NODE CONSTRUCTOR                                 **/
/*********************************************************************************************/
ARGnode::ARGnode( unsigned short int positionOfSourceNode,	unsigned short int positionOfTheMutation , const std::vector<ARGnode>& ARGData ){
	std::vector<unsigned short int>::iterator it;
	std::find(ARGData.at(positionOfSourceNode).mutationPositionList.begin() , ARGData.at(positionOfSourceNode).mutationPositionList.end() , positionOfTheMutation);

	bool mutationPositionPresent = false;
	unsigned short int elementValue;

	for (unsigned short int i = 0 ; i < ARGData.at(positionOfSourceNode).mutationPositionList.size() ;++i){
		elementValue = ARGData.at(positionOfSourceNode).mutationPositionList.at(i);
		if (elementValue != positionOfTheMutation){
			this->mutationPositionList.push_back(elementValue);
		}else{
			mutationPositionPresent = true;
		}
	}

#ifdef VALIDATION
	if (this->mutationPositionList.size() != ARGData.at(positionOfSourceNode).mutationPositionList.size()-1  || !mutationPositionPresent){
		std::cerr << "the mutation position " << positionOfTheMutation << " for the node " << positionOfSourceNode << " is absent." << std::endl;
		std::cerr << "this->mutationPositionList.size() != ARGData.at(positionOfSourceNode).mutationPositionList.size()-1  || !mutationPositionPresent inARGnode::ARGnode ";
		throw  "this->mutationPositionList.size() != ARGData.at(positionOfSourceNode).mutationPositionList.size()-1  || !mutationPositionPresent inARGnode::ARGnode ";
	}
#endif

	this->typeOfNode = ARGNodeType::MUTATIONODE;
	this->ancestralMaterialPositionList = ARGData.at(positionOfSourceNode).ancestralMaterialPositionList;
	this->positionOfDataSource = positionOfSourceNode;
	this->secondDataPoint = positionOfTheMutation;
	this->listOfLeafPositionForTIM = ARGData.at(positionOfSourceNode).listOfLeafPositionForTIM;
}

/*********************************************************************************************/
/**                             COALESCENCE NODE CONSTRUCTOR                                **/
/*********************************************************************************************/
ARGnode::ARGnode( unsigned short int positionOfSourceNode , const std::vector<ARGnode>& ARGData, unsigned short int positionOfTheSecondSouceNode) {

	std::vector<unsigned short int> resultVectorForLeafReach;
#ifdef VALIDATION

#endif

	std::set_intersection( ARGData.at(positionOfSourceNode).ancestralMaterialPositionList.begin(), ARGData.at(positionOfSourceNode).ancestralMaterialPositionList.end(),
			ARGData.at(positionOfTheSecondSouceNode).ancestralMaterialPositionList.begin(),ARGData.at(positionOfTheSecondSouceNode).ancestralMaterialPositionList.end(),
			std::back_inserter(this->ancestralMaterialPositionList));

	std::set_union(ARGData.at(positionOfSourceNode).mutationPositionList.begin(),ARGData.at(positionOfSourceNode).mutationPositionList.end(),
			ARGData.at(positionOfTheSecondSouceNode).mutationPositionList.begin(),ARGData.at(positionOfTheSecondSouceNode).mutationPositionList.end(),
			std::back_inserter(this->mutationPositionList));

	this->typeOfNode = ARGNodeType::COALESCENSENODE;
	this->positionOfDataSource = positionOfSourceNode;
	this->secondDataPoint = positionOfTheSecondSouceNode;

	for (unsigned short int i = 0 ; i < 2*this->sizeOfSequence-2 ; ++i ){ // HIGH COST FOR THIS! CHECK FOR OPTIMIZATION SEMANTICS !!!!!!!!!!
		std::set_union(ARGData.at(positionOfSourceNode).listOfLeafPositionForTIM.at(i).begin() , ARGData.at(positionOfSourceNode).listOfLeafPositionForTIM.at(i).end(),
				ARGData.at(positionOfTheSecondSouceNode).listOfLeafPositionForTIM.at(i).begin() , ARGData.at(positionOfTheSecondSouceNode).listOfLeafPositionForTIM.at(i).end(),
				std::back_inserter(resultVectorForLeafReach));
		this->listOfLeafPositionForTIM.push_back(resultVectorForLeafReach);
		resultVectorForLeafReach.clear();
	}
}

/*********************************************************************************************/
/**                          RECOMBINATION NODE CONSTRUCTOR                                 **/
/*********************************************************************************************/
ARGnode::ARGnode( unsigned short int positionOfSourceNode , unsigned short int positionOfTheRecombination,
		bool trueForleftSideIsAncestralMateriel , const std::vector<ARGnode>& ARGData) {

	std::vector<unsigned short int> tempAncestralMaterialPositionList;
	std::vector<unsigned short int> tempMutationPositionList;

#ifdef VALIDATION
	if ( positionOfTheRecombination >= this->sizeOfSequence ){
		std::cerr << "The position " << positionOfTheRecombination << "is invalid." << std::endl;
		std::cerr << "positionOfTheRecombination >= this->sizeOfSequence" << std::endl;
		throw "positionOfTheRecombination >= this->sizeOfSequence";
	}
	if ( positionOfTheRecombination == 0 ){
		std::cerr << "The position " << positionOfTheRecombination << "is invalid." << std::endl;
		std::cerr << "positionOfTheRecombination == 0 in ARGnode::ARGnode" << std::endl;
		throw "positionOfTheRecombination == 0 in ARGnode::ARGnode";
	}

#endif

	if(trueForleftSideIsAncestralMateriel){
		for (unsigned short int i = 1 ;  i <= positionOfTheRecombination ; ++i){
			tempAncestralMaterialPositionList.push_back(i);

		}
	} else {
		for (unsigned short int i = positionOfTheRecombination+1 ; i <= this->sizeOfSequence ;++i){
			tempAncestralMaterialPositionList.push_back(i);
		}
	}

	std::set_union(tempAncestralMaterialPositionList.begin(), tempAncestralMaterialPositionList.end(),
			ARGData.at(positionOfSourceNode).ancestralMaterialPositionList.begin(), ARGData.at(positionOfSourceNode).ancestralMaterialPositionList.end(),
			std::back_inserter(this->ancestralMaterialPositionList));

	std::set_difference(ARGData.at(positionOfSourceNode).mutationPositionList.begin() , ARGData.at(positionOfSourceNode).mutationPositionList.end(), // Source mutation
			tempAncestralMaterialPositionList.begin(),tempAncestralMaterialPositionList.end(),                                                   // minus new ancestral positions
			std::inserter( this->mutationPositionList , this->mutationPositionList.begin()) );                                                   // in mutations position list

#ifdef VALIDATION
	if ( this->ancestralMaterialPositionList.size() == this->sizeOfSequence){
		std::cerr << "Cant recombinate to form an complete unknown node" << std::endl;
		std::cerr << "this->ancestralMaterialPositionList.size() == this->sizeOfSequence in ARGnode::ARGnode";
		if (this->mutationPositionList.size() != 0){
			std::cerr <<"Serious problem!!!"; // It means that all the material is ancestral but there is still mutations. Normaly this cnat be simulated.
			throw "Serious problem!!! in ARGnode::ARGnode"; // IF this error happens you need to debug IMMEDIALTLY!!
		}
		throw "this->ancestralMaterialPositionList.size() == this->sizeOfSequence in ARGnode::ARGnode";
	}
#endif

	this->typeOfNode = ARGNodeType::RECOMBINATIONODE;
	this->positionOfDataSource = positionOfSourceNode;
	this->secondDataPoint = positionOfTheRecombination;

	std::vector<unsigned short int> emptyVector;
	if(trueForleftSideIsAncestralMateriel){
		for ( unsigned short int i = 0 ; i < 2*this->secondDataPoint-1; ++i){
			this->listOfLeafPositionForTIM.push_back(emptyVector);
		}
		for ( unsigned short int i = 2*this->secondDataPoint-1; i < 2*this->sizeOfSequence-2;++i){
			this->listOfLeafPositionForTIM.push_back(ARGData.at(positionOfSourceNode).listOfLeafPositionForTIM.at(i));
		}
	}else{
		for ( unsigned short int i = 0 ; i < 2*this->secondDataPoint-1; ++i){
			this->listOfLeafPositionForTIM.push_back(ARGData.at(positionOfSourceNode).listOfLeafPositionForTIM.at(i));
		}
		for ( unsigned short int i = 2*this->secondDataPoint-1; i < 2*this->sizeOfSequence-2;++i){
			this->listOfLeafPositionForTIM.push_back(emptyVector);
		}
	}
}

/*********************************************************************************************/
/**                                   DEFAULTED DESTRUCTORS                                **/
/*********************************************************************************************/
/*
ARGnode::~ARGnode() {
	// TODO Auto-generated destructor stub
}
 */


/*********************************************************************************************/
/**                                REDIFINITION OF OUTPUT STREAM                            **/
/*********************************************************************************************/
std::ostream& operator << (std::ostream & out, const ARGnode & aNode){

	if ( aNode.typeOfNode == ARGNodeType::LEAFSAMPLEMODE){
		out << "Leaf node          -" << " with source node origin "<<  aNode.positionOfDataSource+1 << " and with ";
	}else if ( aNode.typeOfNode == ARGNodeType::COALESCENSENODE){
		out << "CoalescenseNode    -" << " with source nodes origin "<<  aNode.positionOfDataSource+1 << " & " << aNode.secondDataPoint+1 << " and with ";
	}else if (aNode.typeOfNode == ARGNodeType::RECOMBINATIONODE){
		out << "Recombination node -" << " with source node origin "<<  aNode.positionOfDataSource+1 << " and with " ;
	}else if ( aNode.typeOfNode == ARGNodeType::MUTATIONODE){
		out << "Mutation node      -" << " with source node origin "<<  aNode.positionOfDataSource+1 << " and with ";
	}

	out << " mutations position at: ";
	for (const auto &i : aNode.mutationPositionList){
		out << i << ' ';
	}
	out << std:: endl;
	out << "   ancestral material position: ";
	for (const auto &i: aNode.ancestralMaterialPositionList){
		out << i << ' ';
	}
	out << std:: endl;
	out << "   Leaf node paths: ";
	for (unsigned short int i = 0 ; i < aNode.listOfLeafPositionForTIM.size() ;++i){
		out << i+1 << ":";
		for(unsigned short int j = 0 ; j < aNode.listOfLeafPositionForTIM.at(i).size() ; ++j){
			out << aNode.listOfLeafPositionForTIM.at(i).at(j)+1;
			if (j < aNode.listOfLeafPositionForTIM.at(i).size() -1 ){
				out<< '-';
			}
		}
		out << "  ";
	}
	return out;
}

