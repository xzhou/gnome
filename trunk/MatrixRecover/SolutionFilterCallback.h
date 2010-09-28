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
	SlnPool sp;
public:
	SolutionFilterCallbackI(IloEnv env) : IloCplex::IncumbentCallbackI(env){};
	void main();	// the call back function

	//the duplicate function to create new call back object
	IloCplex::CallbackI* duplicateCallback() const {
		return (new (getEnv()) SolutionFilterCallbackI(*this));
	}
};

IloCplex::Callback SolutionFilterCallback(IloEnv env);

#endif /* INCUMBENTCALLBACK_H_ */
