//============================================================================
// Name        : SnpSolver.cpp
// Author      : Xiaoyong Zhou
// Version     :
// Copyright   : All Rights Reserved @ zhou@indiana.edu
// Description : MatrixRecover in C, Ansi-style, use IBM cplex to reverse an
//				 binary matrix
//============================================================================

#include <iostream>
#include <fstream>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include "SnpMatrix.h"
#include "AuxFunc.h"
#include "SolutionFilterCallback.h"
#include "ExpConf.h"
#include "Solver.h"
#include "AuxFunc.h"
#include <sstream>
#include "Experiments.h"
int my_debug = 0;

ILOSTLBEGIN

// TODO
//we need a database to remember all unique solutions seen so far

int main(int argc, char *argv[]) {

	string type = "";
	if(argc > 1){
		type = string(argv[1]);
	}

	Experiments exp;

	if(type == "fix"){
		cout<<"fixing and solve"<<endl;
		exp.largeScaleFixRow();
	}
	else if(type == "find"){
		cout<<"find multiple solution"<<endl;
		exp.findMultipleSolution();
	}
	else if(type == "true"){
		//true snps sequence sampling
		//TODO
	}

	return 0;

}
