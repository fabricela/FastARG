/*
 * Fonctions.cpp
 *
 *  Created on: 2015-05-20
 *      Author: Éric Marcotte
 */
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include "Fonctions.h"



/*
 * Check if file opens correctly.
 */
void ouvertureFlux(std::ifstream& fluxFichier,std::string nomFichier){
	fluxFichier.open(nomFichier.c_str());

	if (!fluxFichier.is_open() ){
		std::cerr << "The file:" << nomFichier << "cannot be found!";
		abnormalProgramEnding();
	}
}

/*
 *
 */
SNPSequence getASequenceFromFile(std::ifstream& fluxFichier, int sizeOfSequence){
	int Id = nextInt(fluxFichier);
	bool isCase = nextZeroOrOne(fluxFichier);
	bool isControl = nextZeroOrOne(fluxFichier);

	SNPSequence aSequence(Id,isCase,isControl);

	for (int i = 0 ; i < sizeOfSequence ; ++i){
		if(nextZeroOrOne(fluxFichier)){
			aSequence.mutationPositionList.push_back(i+1);
		}
	}

	/*
	for (unsigned int i = 0 ;i < aSequence.MutationSequence.size(); ++i){
		std::cout << aSequence.MutationSequence[i];
	}
	*/

	return aSequence;
}

/*
 * Exception if an value is invalid !!
 */
class invalidValue: public std::exception
{
	virtual const char* what() const throw(){
		return "One of the value of the file is invalid!";
	}
}invalidValueError;

/*
 *
 */
bool nextZeroOrOne(std::ifstream& fluxFichier){
	char aChar = fluxFichier.get();

	while (aChar == ' ' || aChar == '	' || aChar == '\n'){ // white spaces must be made better... bad programming
		aChar = fluxFichier.get();
	}
	if (aChar == '1'){
		return true;
	} else if (aChar == '0'){
		return false;
	} else {
		throw invalidValueError;
	}
}

int nextInt(std::ifstream&fluxFichier){
	int value=0;
	std::string aStringOfTheValue;
	char aChar = fluxFichier.get();

	while (aChar == ' ' || aChar == '	' || aChar == '\n'){ // white spaces must be made better... bad programming
		aChar = fluxFichier.get();
	}

	while ( aChar >= '0' && aChar <= '9'){
		value += (aChar-48);
		value *= 10;
		aChar = fluxFichier.get();
	}
	value /= 10;

	return value;
}


/*
 * Méthode qui affiche un message signalant la fin du programme.
 */
void NormalProgramEnding(){
	std::cout << "\n\nEXECUTION SUCESSFUL, NO ERRORS DETECTED!" << std::endl;
	exit (EXIT_SUCCESS); // Pour la portabiliter.
}
/*
 * Méthode qui affiche un message signalant la fin du programme.
 */
void abnormalProgramEnding(){
	std::cout << "\n\nERRORS DETECTED!" << std::endl;
	exit(EXIT_FAILURE); //Pour la portabiliter.
}
