/*
 * Solver.h
 *
 *  Created on: Sep 30, 2010
 *      Author: xzhou
 */

#ifndef SOLVER_H_
#define SOLVER_H_

class Solver {

public:
	Solver();
	virtual ~Solver();
	int solve(int **M, int m, int n);	//return the number of solutions
};

#endif /* SOLVER_H_ */
