//============================================================================
// Name        : FastARG.cpp
// Author      : Eric Marcotte
// Version     :
// Copyright   :
// Description : Ancestral recombination graph genereator.
//============================================================================

#define nameOfDataFile     "Example.dat"
#define nameOfparameterFile "Seq_1.4_Ech_3.par"
#define SIZEOFSAMPLE 8
#define SIZEOFSEQUENCE 5

#include <iostream>
#include <algorithm>
#include <limits>
#include <fstream>
#include <random>
#include <vector>
#include <chrono>
#include <bitset>
#include "Fonctions.h"
#include "SNPSequence.h"
#include "ARGraph.h"

typedef std::chrono::high_resolution_clock Clock;

void testExample();
void testExample2();
void testShuffle();
unsigned short int ARGnode::sizeOfSequence = SIZEOFSEQUENCE;
unsigned short int ARGnode::sizeOfSample = SIZEOFSAMPLE;
using std::cout;
using std::endl;

enum class numberOfTIMPresent  {ZERO,ONE};
int main() {

	testExample();
	//testExample2();

	/*
	auto startTime = Clock::now();
	auto endTime = Clock::now();
	std::vector<SNPSequence> data;
	std::ifstream fluxFichier;
	ouvertureFlux(fluxFichier,nameOfDataFile);

	for (int i = 0 ; i < SIZEOFSAMPLE ; ++i ){
		data.push_back(getASequenceFromFile(fluxFichier,SIZEOFSEQUENCE));
	}
	fluxFichier.close();

	startTime = Clock::now();
	ARGraph aARGraph(data,SIZEOFSEQUENCE,SIZEOFSAMPLE);
	endTime = Clock::now();
	std::cout << "ARG initialized in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;
    aARGraph.constructExampleARG();
	aARGraph.printARGSummary();
	aARGraph.generateStatisticData();

	aARGraph.resetARG();
	startTime = Clock::now();
	aARGraph.randomArgGenerator();
	endTime = Clock::now();
	std::cout << "ARG generated in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;
	//aARGraph.generateStatisticData();
	//aARGraph.printARGSummary();
	//aARGraph.printLivesNodes();
	std::system("Pause");
	aARGraph.resetARG();
	 */
	NormalProgramEnding();
	return 0;
}

void testExample(){
	auto startTime = Clock::now();
	auto endTime = Clock::now();
	std::vector<SNPSequence> data;
	std::ifstream fluxFichier;
	ouvertureFlux(fluxFichier,nameOfDataFile);

	for (int i = 0 ; i < SIZEOFSAMPLE ; ++i ){
		data.push_back(getASequenceFromFile(fluxFichier,SIZEOFSEQUENCE));
	}
	fluxFichier.close();

	startTime = Clock::now();
	ARGraph aARGraph(data,SIZEOFSEQUENCE,SIZEOFSAMPLE);
	endTime = Clock::now();
	std::cout << "ARG initialized in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;
	aARGraph.constructExampleARG();
	aARGraph.printARGSummary();
	aARGraph.generateStatisticData();
}

void testExample2(){
	auto startTime = Clock::now();
	auto endTime = Clock::now();
	std::vector<SNPSequence> data;
	std::ifstream fluxFichier;
	ouvertureFlux(fluxFichier,nameOfDataFile);

	for (int i = 0 ; i < SIZEOFSAMPLE ; ++i ){
		data.push_back(getASequenceFromFile(fluxFichier,SIZEOFSEQUENCE));
	}
	fluxFichier.close();

	startTime = Clock::now();
	ARGraph aARGraph(data,SIZEOFSEQUENCE,SIZEOFSAMPLE);
	endTime = Clock::now();
	std::cout << "ARG initialized in: " << std::chrono::duration_cast<std::chrono::milliseconds> (endTime - startTime).count() << " milli seconds" << std::endl;
	aARGraph.randomArgGenerator();
	aARGraph.printARGSummary();
	aARGraph.generateStatisticData();
}

void testShuffle(){
	std::vector<int> vector;
	for (auto i = 1 ; i <= 10 ;++i){
		vector.push_back(i);
	}
	std::random_shuffle(vector.begin(),vector.end());
	for (const auto &i : vector){
		cout<< i << ' ';
	}
}

