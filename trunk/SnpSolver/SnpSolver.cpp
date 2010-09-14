//============================================================================
// Name        : SnpSolver.cpp
// Author      : Xiaoyong Zhou
// Version     :
// Copyright   : All Rights Reserved @ zhou@indiana.edu
// Description : Hello World in C, Ansi-style
//============================================================================

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <readDataToMatrix.h>

ILOSTLBEGIN

int main(void) {
	IloEnv env;
	try {
		IloModel model(env);	//create optimization model
		IloCplex cplex(env);	//create cplex model

		IloObjective obj;
		IloNumVarArray var(env);
		IloRangeArray rng(env);
	}
	catch (IloException& e) {
		cerr << "concert exception caught: " << e << endl;
	}
	env.end();
}
