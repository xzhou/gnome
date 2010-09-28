#include "SolutionFilterCallback.h"
#include <iostream>
#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include "AuxFunc.h"

void SolutionFilterCallbackI::main(){
	//create filter
	//std::cout<< "calling incumbent callback"<<std::endl;
	IloModel model = getModel();
	//cout<<m<<" "<<n<<endl;
	IloEnv env = getEnv();

	int nVar = m*n + m*nchoosek(n, 2);

	IloNumArray val(env, nVar);
	getValues(val, xVar);

	int aSln[m*n];
	for(int i = 0; i < m*n; i ++){
		aSln[i] = val[i];
	}

	if(sp->addToPool(aSln, m ,n) == 0)
	{
		reject();
	}
	//cout<<val<<endl;
}

//export callback
IloCplex::Callback SolutionFilterCallback(IloEnv env, IloNumVarArray x, int m, int n, SlnPool* sp) {
	return (IloCplex::Callback(new (env) SolutionFilterCallbackI(env, x, m, n, sp)));
}
