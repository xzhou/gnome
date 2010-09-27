#include "SolutionFilterCallback.h"
#include <iostream>
#include <stdio.h>


void SolutionFilterCallbackI::main(){
	//create filter
	std::cout<< "calling incumbent callback"<<std::endl;
}

IloCplex::Callback SolutionFilterCallback(IloEnv env) {
	return (IloCplex::Callback(new (env) SolutionFilterCallbackI(env)));
}
