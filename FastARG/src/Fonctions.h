/*
 * Fonctions.h
 *
 *  Created on: 2015-05-18
 *      Author: Éric Marcotte
 */

#ifndef FONCTIONS_H_
#define FONCTIONS_H_
#include <fstream>
#include "SNPSequence.h"

void ouvertureFlux(std::ifstream&,std::string);
SNPSequence getASequenceFromFile(std::ifstream&,int);
int nextInt(std::ifstream&);
bool nextZeroOrOne(std::ifstream&);
void abnormalProgramEnding();
void NormalProgramEnding();

#endif /* FONCTIONS_H_ */
