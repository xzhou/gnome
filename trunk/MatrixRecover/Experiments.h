/*
 * Experiments.h
 * Experiments do different experiments
 *
 *  Created on: Nov 5, 2010
 *      Author: xzhou
 */

#ifndef EXPERIMENTS_H_
#define EXPERIMENTS_H_

#include <fstream>
#include <iostream>
#include "Solver.h"
#include "SnpMatrix.h"
#include "ExpConf.h"

class Experiments {
private:
	ofstream *slnLog;	//solution log
	ofstream *slnCtLog;	//solution count log
	ExpConf conf;
public:
	Experiments();
	virtual ~Experiments();
	void largeScaleFixRow();		//for large scale problem, we verify the difficulty
	void smallScaleFixRow(int m, int n, int repeat);		//
	void findMultipleSolution();
	void exclucde();	//find exclude code
	void testCode();
};

#endif /* EXPERIMENTS_H_ */
