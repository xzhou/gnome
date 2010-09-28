/*
 * IncumbentCallBack.h
 *
 *  Created on: Sep 22, 2010
 *      Author: xzhou
 * implement the incumbent callback functions
 */

#ifndef INCUMBENTCALLBACK_H_
#define INCUMBENTCALLBACK_H_

#include <ilcplex/ilocplex.h>
#include "SlnPool.h"

//export the function

class SolutionFilterCallbackI : public IloCplex::IncumbentCallbackI  {
private:
	//build solution index
	IloNumVarArray xVar;
	int m;
	int n;
	SlnPool *sp;
public:
	SolutionFilterCallbackI	(IloEnv env, IloNumVarArray x, int m, int n, SlnPool* sp) : IloCplex::IncumbentCallbackI(env), xVar(x)
	{
		this->m = m;
		this->n=n;
		this->sp = sp;
	};
	void main();	// the call back function
	//the duplicate function to create new call back object
	IloCplex::CallbackI* duplicateCallback() const {
		return (new (getEnv()) SolutionFilterCallbackI(*this));
	}

};

IloCplex::Callback SolutionFilterCallback(IloEnv env, IloNumVarArray x, int m, int n, SlnPool* sp);

#endif /* INCUMBENTCALLBACK_H_ */
