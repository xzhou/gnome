/*
 * Solver.h
 *
 *  Created on: Sep 30, 2010
 *      Author: xzhou
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <fstream>
#include "SnpMatrix.h"
using namespace std;


ILOSTLBEGIN

class Solver {
private:
	ofstream logs;
public:
	Solver(string logFileName = "solver.log");
	virtual ~Solver();
	int solve(SnpMatrix &M);	//return the number of solutions, use incumbent call back
	int solveAndFilter(SnpMatrix &M);	//solve, store all solution and filter
};

#endif /* SOLVER_H_ */
