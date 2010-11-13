/*
 * Solver.h
 *
 *  Created on: Sep 30, 2010
 *      Author: xzhou
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include <fstream>
#include "SlnPool.h"
#include "SnpMatrix.h"
using namespace std;

#define POP_LIM 250

ILOSTLBEGIN

class Solver {
private:
	ofstream logs;
public:
	Solver(string logFileName = "solver.log");
	virtual ~Solver();
	SlnPool* solveWithCallback(SnpMatrix &M);	//return the number of solutions, use incumbent call back
	SlnPool* solveAndFilter(SnpMatrix &M);	//solve, store all solution and filter
	SlnPool* fixAndSolve(SnpMatrix &M, const int rowIndex);	//fix a row in M
	SlnPool* fixAndSolve(SnpMatrix &M, const int* fixVar);	//fix arbitrary

};

#endif /* SOLVER_H_ */
