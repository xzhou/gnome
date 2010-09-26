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

class SolutionFilterCallback : public IloCplex::IncumbentCallbackI {
public:
	void main();	// the call back function

	//the duplicate function to create new call back object
	IloCplex::CallbackI* duplicateCallback() const {
		return (new (getEnv()) MyCallbackI(*this));
	}
}


#endif /* INCUMBENTCALLBACK_H_ */
