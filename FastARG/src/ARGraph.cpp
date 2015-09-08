/*
 * ARGraph.cpp
 *
 *  Created on: 2015-09-05
 *      Author: Éric Marcotte
 */

#include "ARGraph.h"
#include "ARGnode.h"
#include <algorithm>
#include <limits>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#define VALIDATION
//#define EXECUTIONTRACE

/*********************************************************************************************/
/**                     ANCESTRAL RECOMBINATION GRAPH CONSTRUCTOR                           **/
/*********************************************************************************************/
ARGraph::ARGraph(const std::vector<SNPSequence>& dataSequences,unsigned short int sizeOfSequence,unsigned short int sizeOfSample) {
	this->sizeOfSample = sizeOfSample;
	this->sizeOfSequence = sizeOfSequence;
	this->dataSequences = dataSequences;

#ifdef EXECUTIONTRACE
	std::cout << "dataSeq:" << std::endl;
	for(unsigned int i = 0 ; i < this->dataSequences.size();++i ){
		std::cout << i << ":";
		for (unsigned short int j = 0 ; j< this->dataSequences.at(i).mutationPositionList.size() ;++j){
			std::cout << this->dataSequences.at(i).mutationPositionList.at(j) << ' ';
		}
		std::cout << std::endl;
	}
#endif

#ifdef VALIDATION
	if ( this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() ){
		std::cerr << "this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() in ARGraph::ARGraph";
		throw "this->sizeOfSample >= std::numeric_limits<unsigned short int>::max() in ARGraph::ARGraph";
	}
	if (this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() ){
		std::cerr << "this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() in ARGraph::ARGraph";
		throw "this->sizeOfSequence >= std::numeric_limits<unsigned short int>::max() in ARGraph::ARGraph";
	}
#endif

	unsigned short int i = 0;
	for ( const auto &aData : this->dataSequences){
		ARGnode aLeafnode( aData.mutationPositionList , i);
		this->ARGnodes.push_back(aLeafnode);
		this->positionOfLiveNodes.push_back(i);
		if(aData.asTIM){
			this->leafPositionWithTheTIM.push_back(i);
		}

		if(  (i%2 == 0) ){
			if (this->dataSequences.at(i).isACase !=  this->dataSequences.at(i+1).isACase ){
				std::cerr << "inconsitent data input at :" << i << " & " << i+1 << " leafs nodes." << std::endl;
				std::cerr << "this->dataSequences.at(i).isACase !=  this->dataSequences.at(i+1).isACase in ARGraph::ARGraph";
				throw "this->dataSequences.at(i).isACase !=  this->dataSequences.at(i+1).isACase in ARGraph::ARGraph";
			}
			if (aData.isACase){
				this->sampleIDThatAreCases.push_back(aData.sampleId );
			}
		}
		++i;
	}

#ifdef EXECUTIONTRACE

	for (const auto &i : this->ARGnodes ){
		std::cout << "test :" << i;
	}
#endif
}



/*********************************************************************************************/
/**                                                           **/
/*********************************************************************************************/
void ARGraph::generateACoalescentEvent(unsigned short int firstSourceNode , unsigned short int secondSourceNode) {

	std::vector<unsigned short int> resultOfSetSymmetricDifference;
	std::vector<unsigned short int> resultOfSetUnion;
	std::vector<unsigned short int> result;

#ifdef VALIDATION
	if (firstSourceNode == secondSourceNode){
		std::cerr << "The node sources cant be the same in a coalescent event." << std::endl;
		std::cerr << "The node sources cant be the same in a coalescent event." ;
		throw "The node sources cant be the same in a coalescent event.";
	}
#endif

	if (firstSourceNode > secondSourceNode){
		unsigned short int temp = firstSourceNode;
		firstSourceNode = secondSourceNode;
		secondSourceNode = temp;
	}
	std::vector<unsigned short int>::iterator it1;
	it1 = std::find(this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end() , firstSourceNode);
	if ( it1 == this->positionOfLiveNodes.end() ){
		std::cerr << "Position " << firstSourceNode+1 << " is not a live node." << std::endl;
		std::cerr << "it == this->positionOfLiveNodes.end() in ARGraph::generateACoalescentEvent";
		throw "it == this->positionOfLiveNodes.end() in ARGraph::generateACoalescentEvent";
	}
	std::vector<unsigned short int>::iterator it2;
	it2 = std::find(this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end() , secondSourceNode);
	if ( it2 == this->positionOfLiveNodes.end() ){
		std::cerr << "Position " << secondSourceNode+1 << " is not a live node." << std::endl;
		std::cerr << "it == this->positionOfLiveNodes.end() in ARGraph::generateACoalescentEvent";
		throw "it == this->positionOfLiveNodes.end() in ARGraph::generateACoalescentEvent";
	}

#ifdef VALIDATION
	// Check if this computation is necessary or suffisant.
	std::set_symmetric_difference(this->ARGnodes.at(firstSourceNode).mutationPositionList.begin() , this->ARGnodes.at(firstSourceNode).mutationPositionList.end(),
			this->ARGnodes.at(secondSourceNode).mutationPositionList.begin() , this->ARGnodes.at(secondSourceNode).mutationPositionList.end(),
			std::back_inserter(resultOfSetSymmetricDifference));


	std::set_union( this->ARGnodes.at(firstSourceNode).ancestralMaterialPositionList.begin() , this->ARGnodes.at(firstSourceNode).ancestralMaterialPositionList.end(),
			this->ARGnodes.at(secondSourceNode).ancestralMaterialPositionList.begin()  , this->ARGnodes.at(secondSourceNode).ancestralMaterialPositionList.end(),
			std::back_inserter(resultOfSetUnion));

	std::set_difference(resultOfSetSymmetricDifference.begin(), resultOfSetSymmetricDifference.end(),
			resultOfSetUnion.begin(),resultOfSetUnion.end(),
			std::inserter( result,result.begin()));

	if( result.size() != 0 ){
		std::cerr << "result.size() != 0 in ARG::generateACoalescentEvent";
		throw "result.size() != 0 in ARG::generateACoalescentEvent";
	}

#endif

	ARGnode aCoalescenceNode( firstSourceNode, this->ARGnodes , secondSourceNode );
	this->ARGnodes.push_back(aCoalescenceNode);
	(*it1) = this->ARGnodes.size()-1;
	this->positionOfLiveNodes.erase(it2);
}

/*********************************************************************************************/
/**                                                           **/
/*********************************************************************************************/
void ARGraph::generateARecombinationEvent(unsigned short int positionOfTheRecombination , unsigned short int sourceNode) {

	std::vector<unsigned short int>::iterator it;
	it = std::find(this->positionOfLiveNodes.begin() , this->positionOfLiveNodes.end(),sourceNode  );

#ifdef VALIDATION
	if (it == this->positionOfLiveNodes.end() ){
		std::cerr << "Position " << sourceNode+1 << " is not a live node." << std::endl;
		std::cerr << "it == this->positionOfLiveNodes.end() inARGraph::generateARecombinationEvent";
		throw "it == this->positionOfLiveNodes.end() inARGraph::generateARecombinationEvent";
	}

#endif

	ARGnode leftNode (sourceNode,positionOfTheRecombination,true,this->ARGnodes);
	ARGnode rightNode (sourceNode,positionOfTheRecombination,false,this->ARGnodes);
	this->ARGnodes.push_back(leftNode);
	(*it) = this->ARGnodes.size()-1;
	this->ARGnodes.push_back(rightNode);
	this->positionOfLiveNodes.push_back(this->ARGnodes.size()-1);
}

void ARGraph::generetaAMutationEvent(unsigned short int positionOfTheMutation, unsigned short int sourceNode) {
	std::vector<unsigned short int>::iterator it;
	it = std::find(this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end(),sourceNode);

#ifdef VALIDATION
	if (it == this->positionOfLiveNodes.end()){
		std::cerr << sourceNode << " is not a live node." << std::endl;
		std::cerr << sourceNode << "it == this->positionOfLiveNodes.end() in ARGraph::generetaAMutationEvent." << std::endl;
		throw  "it == this->positionOfLiveNodes.end() in ARGraph::generetaAMutationEvent.";
	}
	if ( positionOfTheMutation > this->sizeOfSequence){
		std::cerr << "The position " << positionOfTheMutation << " is invalid." << std::endl;
		std::cerr << "positionOfTheMutation >= this->sizeOfSequence in  ARGraph::generetaAMutationEvent" << std::endl;
		throw "positionOfTheMutation >= this->sizeOfSequence in  ARGraph::generetaAMutationEvent";
	}
#endif
	ARGnode aMutationNode( sourceNode, positionOfTheMutation,this->ARGnodes );
	this->ARGnodes.push_back(aMutationNode);
	(*it) =  this->ARGnodes.size()-1;

}

/*********************************************************************************************/
/**                                                          **/
/*********************************************************************************************/
void ARGraph::printLivesNodes() {
	std::cout << this->positionOfLiveNodes.size() << " nodes lives with thoses position in ARGnodes vector: ";
	for (unsigned short int i = 0 ; i < this->positionOfLiveNodes.size(); ++i){
		std::cout << this->positionOfLiveNodes.at(i)+1 << ' ';
	}
	std::cout << std::endl;

	for(const auto &i : this->positionOfLiveNodes){
		std::cout << '(' << i+1 << ") "<< this->ARGnodes.at(i) << std::endl;
	}
	std::cout << std::endl;
}

// Destructors
/*
ARGraph::~ARGraph() {
	// TODO Auto-generated destructor stub
}
 */

/*********************************************************************************************/
/**  THIS FUNCTION CONSTUCT A ARG WITH A SPECIFIC SET OF DATA FOR TESTING PURPOSES          **/
/**  HERE IS THE DATA YOU NEED TO INPUT TO GET THE MRCA:                                    **/
/*********************************************************************************************/
void ARGraph::constructExampleARG() {
	auto startTime = Clock::now();
	auto endTime = Clock::now();

	//std::random_shuffle(this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end()); // Should not change anything.
	printLivesNodes();
	generateACoalescentEvent(6,7);
	printLivesNodes();
	generateARecombinationEvent(2,5);
	printLivesNodes();
	generateARecombinationEvent(1,4);
	printLivesNodes();
	generateARecombinationEvent(2,2);
	printLivesNodes();
	generateACoalescentEvent(1,14);
	printLivesNodes();
	generetaAMutationEvent(2 ,15);
	printLivesNodes();
	generateACoalescentEvent(3,13);
	printLivesNodes();
	generateACoalescentEvent(12,17);
	printLivesNodes();
	generateACoalescentEvent(10,11);
	printLivesNodes();
	generetaAMutationEvent(5,8);
	printLivesNodes();
	generateACoalescentEvent(9,20);
	printLivesNodes();
	generateARecombinationEvent(1,18);
	printLivesNodes();
	generateACoalescentEvent(23,16);
	printLivesNodes();
	generateACoalescentEvent(0,24);
	printLivesNodes();
	generetaAMutationEvent(1,25);
	printLivesNodes();
	generetaAMutationEvent(4 ,19);
	printLivesNodes();
	generateACoalescentEvent(22,27);
	printLivesNodes();
	generateACoalescentEvent(28,21);
	printLivesNodes();
	generetaAMutationEvent(3 ,29);
	printLivesNodes();
	generateACoalescentEvent(26,30);
	printLivesNodes();
	this->positionOfLiveNodes.clear();
	endTime = Clock::now();
	std::cout << "ARG generated in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;
}

void ARGraph::generateStatisticData() {
	auto startTime = Clock::now();
	auto endTime = Clock::now();
	std::vector<unsigned short int> intersection;
	int numberOfmatch = 0;

	for (unsigned short int i = 0 ; i < 2*this->sizeOfSequence-2;++i){
		for(unsigned short int j = 0; j < this->ARGnodes.size()-1 ;++j ){
			std::set_intersection(this->ARGnodes.at(j).listOfLeafPositionForTIM.at(i).begin(), this->ARGnodes.at(j).listOfLeafPositionForTIM.at(i).end(),
					this->leafPositionWithTheTIM.begin(),this->leafPositionWithTheTIM.end(),
					std::back_inserter(intersection));
			numberOfmatch += intersection.size()*intersection.size();
			intersection.clear();
		}
		this->statistics.push_back(numberOfmatch);
		numberOfmatch = 0;
	}
	std::cout << "Numbers Of match by position:" << std::endl;
	unsigned short int position = 1;
	for (const auto &i : this->statistics){
		std::cout << position << ": " << i << std::endl;
		++position;
	}
	endTime = Clock::now();
	std::cout << "Statistics generated in: "
			<< std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count()
			<< " milliseconds" << std::endl;
}

bool ARGraph::generateARandomCoalescentEvent() {
	std::vector<unsigned short int> resultOfSetSymmetricDifference;
	std::vector<unsigned short int> resultOfSetUnion;
	std::vector<unsigned short int> result;

	for (unsigned short int i = 0 ; i < this->positionOfLiveNodes.size()-1 ;++i){
		for (unsigned short int j = i+1 ; j < this->positionOfLiveNodes.size() ;++j){
			if( this->ARGnodes.at(this->positionOfLiveNodes.at(i)).mutationPositionList.size() ==
					this->ARGnodes.at(this->positionOfLiveNodes.at(j)).mutationPositionList.size()) {

				std::set_symmetric_difference(
						this->ARGnodes.at(this->positionOfLiveNodes.at(i)).mutationPositionList.begin(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(i)).mutationPositionList.end(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(j)).mutationPositionList.begin(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(j)).mutationPositionList.end(),
						std::back_inserter(resultOfSetSymmetricDifference));


				std::set_union(
						this->ARGnodes.at(this->positionOfLiveNodes.at(i)).ancestralMaterialPositionList.begin(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(i)).ancestralMaterialPositionList.end(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(j)).ancestralMaterialPositionList.begin(),
						this->ARGnodes.at(this->positionOfLiveNodes.at(j)).ancestralMaterialPositionList.end(),
						std::back_inserter(resultOfSetUnion));

				std::set_difference(
						resultOfSetSymmetricDifference.begin(),
						resultOfSetSymmetricDifference.end(),
						resultOfSetUnion.begin(),
						resultOfSetUnion.end(),
						std::inserter( result,result.begin()));

#ifdef EXECUTIONTRACE
				std::cout << "result for i=" << i << " j=" << j << " : ";
				for (const auto &k : result){
					std::cout << k << ' ';
				}
				std::cout << std::endl;
#endif
				if ( result.size() == 0 ){
					std::cout << "Found a coalescense between node " << this->positionOfLiveNodes.at(i)+1 << " and " << this->positionOfLiveNodes.at(j)+1 << std::endl;
					generateACoalescentEvent(this->positionOfLiveNodes.at(i),this->positionOfLiveNodes.at(j));
					return true;
				}
				resultOfSetSymmetricDifference.clear();
				resultOfSetUnion.clear();
				result.clear();
			}
		}
	}
	return false;
}

void ARGraph::randomArgGenerator() {

	//check if empty
	auto startTime = Clock::now();
	auto endTime = Clock::now();
/*
	while ( this->positionOfLiveNodes.size()  != 1){
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomMutationEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomRecombinationEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::cout << " Size of lives nodes vector: " <<  this->positionOfLiveNodes.size() << std::endl;
		std::system("pause");
	}
	*/

	while ( ( this->positionOfLiveNodes.size() != 1) ||
			(this->ARGnodes.at(this->positionOfLiveNodes.at(0)).mutationPositionList.size() != 0) ){

		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomMutationEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomRecombinationEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomCoalescentEvent();
		std::random_shuffle( this->positionOfLiveNodes.begin(),this->positionOfLiveNodes.end());
		generateARandomMutationEvent();
	}
	endTime = Clock::now();
	std::cout << "ARG generated in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;

}

void ARGraph::resetARG() {
	this->statistics.clear();
	this->positionOfLiveNodes.clear();
	for (unsigned short int i = 0 ;i < this->sizeOfSample; ++i){
		this->positionOfLiveNodes.push_back(i);
	}
	this->ARGnodes.erase(this->ARGnodes.begin()+8,this->ARGnodes.end());
}

bool ARGraph::generateARandomMutationEvent() {

	std::mt19937 eng((std::random_device())());
	std::uint32_t seed_val;
	eng.seed(seed_val);
	std::uniform_int_distribution<> dist(0, (this->positionOfLiveNodes.size()-1));

	/*
	while ( aRandomInt != this->positionOfLiveNodes.size()-1 ){
		std::cout << aRandomInt << " ";
		aRandomInt = dist(eng);
	}
	std::cout << aRandomInt;
	 */

	/*
	int nombreDe1 =0;
	for (int i = 0 ; i < 100000;++i){
		int aRandomInt = dist(eng);
		std::cout << aRandomInt;
		if (aRandomInt == 1){
			++nombreDe1;
		}
	}
	std::cout << std::endl << "Pourcentage de 1 = " << (float) nombreDe1/100000 << "%.";
	 */

	int randomNode = dist(eng);
	int randomMutation;
	int numberOfMutationForTheRandomNode = this->ARGnodes.at(this->positionOfLiveNodes.at(randomNode)).mutationPositionList.size();
	int mutation;
	if ( numberOfMutationForTheRandomNode != 0){

		std::uniform_int_distribution<> dist2(1,  numberOfMutationForTheRandomNode );
		randomMutation = dist2(eng);
		mutation =  this->ARGnodes.at(this->positionOfLiveNodes.at(randomNode)).mutationPositionList.at(randomMutation-1);

		std::cout << "Mutation for node " << this->positionOfLiveNodes.at(randomNode)+1 << " for mutation position: "
				<< mutation  << "."<< std::endl;
		generetaAMutationEvent( mutation , this->positionOfLiveNodes.at(randomNode));
		return true;
	}else{
		return false;
	}
}

bool ARGraph::generateARandomRecombinationEvent() {
	std::mt19937 eng((std::random_device())());
	std::uint32_t seed_val;
	eng.seed(seed_val);
	std::uniform_int_distribution<> dist(0, (this->positionOfLiveNodes.size()-1));

	int randomNode = this->positionOfLiveNodes.at(dist(eng));
	if (this->ARGnodes.at(randomNode).mutationPositionList.size() < 2 ){
		return false;
	}else{
		std::uniform_int_distribution<> dist2( 1, this->sizeOfSequence-1);
		try{
			int temp = dist2(eng);
			std::cout << "Recombination event with point: " << temp << " at node:" << randomNode+1 << '.' << std::endl;
			generateARecombinationEvent(dist2(eng),randomNode);

			return true;
		}catch(...){
			return false;
		}

	}
}

void ARGraph::printARGSummary() {
	std::cout << "This ARG has " << this->positionOfLiveNodes.size() << " lives nodes." << std::endl;
	std::cout << " List of leaf position with the TIM: ";
	for (auto const &i : this->leafPositionWithTheTIM ){
		std::cout << i+1 << ' ';
	}
	std::cout << std::endl;
	std::cout << " List of cases in sample IDs: ";
	for (auto const &i : this->sampleIDThatAreCases ){
		std::cout << i << ' ';
	}
	std::cout << std::endl;



}
