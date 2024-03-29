/*
 * Solver.cpp
 *
 *  Created on: Sep 30, 2010
 *      Author: xzhou
 */

#include <fstream>
#include <ilcplex/ilocplex.h>
#include <string>
#include "Solver.h"
#include "AuxFunc.h"
#include "SlnPool.h"
#include "SnpMatrix.h"
#include "SolutionFilterCallback.h"

Solver::Solver(string logFileName) :
logs(logFileName.c_str()) {
	// TODO Auto-generated constructor stub
}

Solver::~Solver() {
	// TODO Auto-generated destructor stub
}


//this function return the number solutions for M
SlnPool * Solver::solveWithCallback(SnpMatrix &M){
	//SnpMatrix M = SnpMatrix(M_in, m, n);
	//create model, calculate the number of variables
	int m = M.nInd;
	int n = M.nSnp;
	int nPair = nchoosek(n, 2);
	int nAuxVar = m*nPair;
	int nVar = m*n + nAuxVar;

	int numSolution = -1;

	SlnPool *sp = new SlnPool();

	IloEnv env;
	ofstream c1("c1.log");
	env.setOut(c1);
	try {
		IloModel model(env);	//create optimization model

		//initialize the model
		IloNumArray obj(env);
		IloNumVarArray x(env);

		//initialize the obj function
		for (int i = 0; i < nVar; i ++){
			obj.add(0);
			x.add(IloNumVar(env, 0, 1, ILOINT));
		}
		if (obj.getSize() != nVar) {
			cerr<<"initializing the target"<<endl;
		}
		if (x.getSize() != obj.getSize()){
			cerr<<"X:"<<x.getSize()<<" obj: "<<obj.getSize()<<endl;
		}

		//create objective as 0
		model.add(IloMaximize(env, 0));

		IloRangeArray con(env);
		//add p constraints
		for (int j = 0; j < n; j ++){
			IloExpr pcons(env);
			for(int i = 0; i < m;  i ++){
				pcons += x[j*m + i] * 1;
			}
			con.add(M.P[j] <= pcons <= M.P[j]);
		}
		cout<<"p constraints added"<<endl;
		//add R constraints, e.g. c00 constraints
		int base = m*n;
		for(int i = 0; i < nPair; i ++){
			int begin = m*i;
			int end = m*(i+1)-1;

			IloExpr rcons(env);
			for(int j = begin; j <= end; j ++){
				rcons += x[base + j] * 1;
			}
			con.add(M.LR[i] <= rcons <= M.LR[i]);
		}
		//cout<<"r constraints added"<<endl;

		//add x_ij ~ axu_ijk constraints
		for (int j = 0; j < n-1; j ++) {
			for (int k = j+1; k < n; k ++ ) {
				for (int i = 0; i < m; i ++ ) {
					int d1 = j*m + i;
					int d2 = k*m + i;
					int dd_relative = calcdd(i, j, k, m, n);
					int dd = dd_relative + m*n - 1;
					//cout<<d1<<" "<<d2<<" "<<dd_relative<<"/"<<dd<<endl;
					IloExpr auxConsTop(env);
					auxConsTop = - x[d1] - x[d2] + 2 * x[dd];
					con.add(-INFINITY <= auxConsTop <= 0);

					IloExpr auxConsBtn(env);
					auxConsBtn =   x[d1] + x[d2] - 2 * x[dd];
					con.add(-INFINITY <= auxConsBtn <= 1);
				}
			}
		}
		//cout<<"pairwise constraint added"<<endl;
		model.add(con);
		//export the model for verification
		IloCplex cplex(model);
		cplex.exportModel("r.lp");

		//set parameters before populate solutions
		cplex.setParam(IloCplex::SolnPoolIntensity, 4);	//enum all solutions
		cplex.setParam(IloCplex::SolnPoolGap, 0);
		IloCplex::Callback mycallback = cplex.use(SolutionFilterCallback(env, x, m, n, sp));
		cplex.populate();

		//check the solutions
		numSolution = cplex.getSolnPoolNsolns();
		env.out()<<"the solution pool contains "<<numSolution<<" solutions"<<endl;
		mycallback.end();
		env.end();
		//sp->printPool(cout);
	}
	catch (IloException& e) {
		cerr << "concert exception caught: " << e << endl;
	}
	c1.close();
	env.end();
	return sp;
}

//find experiments
SlnPool * Solver::solveAndFilter(SnpMatrix &M){
	//SnpMatrix M = SnpMatrix(M_in, m, n);
	//create model, calculate the number of variables
	int m = M.nInd;
	int n = M.nSnp;
	int nPair = nchoosek(n, 2);
	int nAuxVar = m*nPair;
	int nVar = m*n + nAuxVar;

	int numSolution = -1;
	int nRealSln = 0;
	SlnPool *sp = new SlnPool();

	IloEnv env;
	ofstream c2("c2.log");
	env.setOut(c2);
	try {
		IloModel model(env);	//create optimization model

		//initialize the model
		IloNumArray obj(env);
		IloNumVarArray x(env);

		//initialize the obj function
		for (int i = 0; i < nVar; i ++){
			obj.add(0);
			x.add(IloNumVar(env, 0, 1, ILOINT));
		}
		if (obj.getSize() != nVar) {
			cerr<<"initializing the target"<<endl;
		}
		if (x.getSize() != obj.getSize()){
			cerr<<"X:"<<x.getSize()<<" obj: "<<obj.getSize()<<endl;
		}

		//create objective as 0
		model.add(IloMaximize(env, 0));

		IloRangeArray con(env);
		//add p constraints
		for (int j = 0; j < n; j ++){
			IloExpr pcons(env);
			for(int i = 0; i < m;  i ++){
				pcons += x[j*m + i] * 1;
			}
			con.add(M.P[j] <= pcons <= M.P[j]);
		}
		//cout<<"p constraints added"<<endl;
		//add R constraints, e.g. c00 constraints
		int base = m*n;
		for(int i = 0; i < nPair; i ++){
			int begin = m*i;
			int end = m*(i+1)-1;

			IloExpr rcons(env);
			for(int j = begin; j <= end; j ++){
				rcons += x[base + j] * 1;
			}
			con.add(M.LR[i] <= rcons <= M.LR[i]);
		}
		//cout<<"r constraints added"<<endl;

		//add x_ij ~ axu_ijk constraints
		for (int j = 0; j < n-1; j ++) {
			for (int k = j+1; k < n; k ++ ) {
				for (int i = 0; i < m; i ++ ) {
					int d1 = j*m + i;
					int d2 = k*m + i;
					int dd_relative = calcdd(i, j, k, m, n);
					int dd = dd_relative + m*n - 1;
					//cout<<d1<<" "<<d2<<" "<<dd_relative<<"/"<<dd<<endl;
					IloExpr auxConsTop(env);
					auxConsTop = - x[d1] - x[d2] + 2 * x[dd];
					con.add(-INFINITY <= auxConsTop <= 0);

					IloExpr auxConsBtn(env);
					auxConsBtn =   x[d1] + x[d2] - 2 * x[dd];
					con.add(-INFINITY <= auxConsBtn <= 1);
				}
			}
		}

		model.add(con);
		//export the model for verification
		IloCplex cplex(model);
		cplex.exportModel("r.lp");

		//set parameters before populate solutions
		cplex.setParam(IloCplex::SolnPoolIntensity, 4);	//enum all solutions
		cplex.setParam(IloCplex::SolnPoolAGap, 0.0);
		cplex.setParam(IloCplex::PopulateLim, POP_LIM);	//the top of solutions
//		cplex.setParam(IloCplex::PopulateSolLim, );
//		cplex.setParam(IloCplex::SolnPoolCapacity, 10000000000);
		//IloCplex::Callback mycallback = cplex.use(SolutionFilterCallback(env, x, m, n, sp));
		cplex.populate();

		//check the solutions
		numSolution = cplex.getSolnPoolNsolns();
		//TODO access solutions and remove duplicate solutions
		for(int si = 0; si < numSolution; si++){
			IloNumArray vals(env);
			cplex.getValues(vals, x, si);
			if(sp->addToPool(vals, m, n)){
				nRealSln ++;
			}
		}

		//sp->printPool(cout);

		env.out()<<"the solution pool contains "<<numSolution<<" solutions"<<endl;
		env.end();
	}
	catch (IloException& e) {
		cerr << "concert exception caught: " << e << endl;
	}
	c2.close();
	env.end();
	return sp;
}

SlnPool* Solver::fixAndSolve(SnpMatrix &M, const int fixRow){
	int *fixVar = M.M[fixRow];

	return fixAndSolve(M, fixVar);
}

/*
 *	This method try to prove that given a solution,
 *  the attacker can not fix one of them using fix and solve method
 *	This
 */
SlnPool* Solver::fixAndSolve(SnpMatrix &M, const int* fixVar){
	int m = M.nInd;
	int n = M.nSnp;

	//int *fixVar = M.M[fixRow];

	int nPair = nchoosek(n, 2);		//number of pairwise constraints
	int nAuxVar = m*nPair;			//auxiliary variable
	int nFixVar = m*n;				//fixing one row we need this
	int nVar = m*n + nFixVar + nAuxVar;

	int numSolution = -1;
	//int nRealSln = 0;
	SlnPool *sp = new SlnPool();

	IloEnv env;
	ofstream c2("c2.log");
	env.setOut(c2);

	try{
		IloModel model(env);

		//initialize the obj function
		IloNumArray obj(env);
		IloNumVarArray x(env);

		//initialize the obj function
		for (int i = 0; i < nVar; i ++){
			obj.add(0);
			x.add(IloNumVar(env, 0, 1, ILOINT));
		}
		if (obj.getSize() != nVar) {
			cerr<<"initializing the target"<<endl;
		}
		if (x.getSize() != obj.getSize()){
			cerr<<"X:"<<x.getSize()<<" obj: "<<obj.getSize()<<endl;
		}

		//create objective as 0
		model.add(IloMaximize(env, 0));

		IloRangeArray con(env);
		//add p constraints
		for (int j = 0; j < n; j ++){
			IloExpr pcons(env);
			for(int i = 0; i < m;  i ++){
				pcons += x[j*m + i] * 1;
			}
			con.add(M.P[j] <= pcons <= M.P[j]);
		}
		//cout<<"p constraints added"<<endl;
		//add R constraints, e.g. c00 constraints
		int base = m*n;
		for(int i = 0; i < nPair; i ++){
			int begin = m*i;
			int end = m*(i+1)-1;

			IloExpr rcons(env);
			for(int j = begin; j <= end; j ++){
				rcons += x[base + j] * 1;
			}
			con.add(M.LR[i] <= rcons <= M.LR[i]);
		}
		//cout<<"r constraints added"<<endl;

		//add x_ij ~ axu_ijk constraints
		for (int j = 0; j < n-1; j ++) {
			for (int k = j+1; k < n; k ++ ) {
				for (int i = 0; i < m; i ++ ) {
					int d1 = j*m + i;
					int d2 = k*m + i;
					int dd_relative = calcdd(i, j, k, m, n);
					int dd = dd_relative + m*n - 1;
					//cout<<d1<<" "<<d2<<" "<<dd_relative<<"/"<<dd<<endl;
					IloExpr auxConsTop(env);
					auxConsTop = - x[d1] - x[d2] + 2 * x[dd];
					con.add(-INFINITY <= auxConsTop <= 0);

					IloExpr auxConsBtn(env);
					auxConsBtn =   x[d1] + x[d2] - 2 * x[dd];
					con.add(-INFINITY <= auxConsBtn <= 1);
				}
			}
		}
		//adding fixing constraints
		int fixVarBeginIndex = m*n + nAuxVar;	//beginning of aux variable
		for(int i = 0; i < m; i ++){		//for each row
			for(int j = 0; j < n; j ++){	//for each col
				int x_ij = j*m+i;	//index variable
				int fixVarIdx = fixVarBeginIndex + j*m+i;

				IloExpr fixVarCons1(env);
				fixVarCons1 = x[fixVarIdx] - x[x_ij];
				con.add(-INFINITY <= fixVarCons1 <= fixVar[j]);

				IloExpr fixVarCons2(env);
				fixVarCons2 = x[fixVarIdx] + x[x_ij];
				con.add(-INFINITY <= fixVarCons2 <= 2 - fixVar[j]);

				IloExpr fixVarCons3(env);
				fixVarCons3 = x[fixVarIdx] + x[x_ij];
				con.add(fixVar[j] <= fixVarCons3 <= INFINITY);

				IloExpr fixVarCons4(env);
				fixVarCons4 = x[fixVarIdx] - x[x_ij];
				con.add(-fixVar[j] <= fixVarCons4 <= INFINITY);
			}
		}

		//constraints for fixing variable
		for(int i = 0; i < m; i ++){
			IloExpr fcons(env);
			for(int j = 0; j < n; j ++){
				int idx = fixVarBeginIndex + j*m+i;
				fcons += x[idx]*1;
			}
			con.add(1 <= fcons <= n);	//for each row, at least one bit different.
		}
		model.add(con);

		IloCplex cplex(model);
		cplex.exportModel("fixSolve.lp");

		//set parameters before populate solutions
		cplex.setParam(IloCplex::SolnPoolIntensity, 4);	//enum all solutions
		cplex.setParam(IloCplex::SolnPoolAGap, 0.0);
		cplex.setParam(IloCplex::PopulateLim, POP_LIM);	//the top of solutions
//		cplex.setParam(IloCplex::PopulateSolLim, );
//		cplex.setParam(IloCplex::SolnPoolCapacity, 10000000000);
		cplex.populate();

		//check the solutions
		numSolution = cplex.getSolnPoolNsolns();
		int nRealSln = 0;
		//TODO access solutions and remove duplicate solutions
		for(int si = 0; si < numSolution; si++){
			IloNumArray vals(env);
			cplex.getValues(vals, x, si);
			if(sp->addToPool(vals, m, n)){
				nRealSln ++;
			}
		}

		//sp->printPool(cout);

		env.out()<<"the solution pool contains "<<numSolution<<" solutions"<<endl;
		env.end();
	}
	catch (IloException &e)
	{
		cerr<<"concert exception caught: "<< e << endl;
	}

	env.end();
	return sp;
}

//This function
SlnPool* Solver::excludeAll(SnpMatrix &M){
	int m = M.nInd;
	int n = M.nSnp;

	int nPair = nchoosek(n, 2);		//number of pairwise constraints
	int nAuxVar = m*nPair;			//auxiliary variable
	int nFixVar = m*n * m;				//fix all row
	int nVar = m*n + nFixVar + nAuxVar;

	int numSolution = -1;
	//int nRealSln = 0;
	SlnPool *sp = new SlnPool();

	IloEnv env;
	ofstream c2("exclude.log");
	env.setOut(c2);

	try{
		IloModel model(env);

		//initialize the obj function
		IloNumArray obj(env);
		IloNumVarArray x(env);

		//initialize the obj function
		for (int i = 0; i < nVar; i ++){
			obj.add(0);
			x.add(IloNumVar(env, 0, 1, ILOINT));
		}
		if (obj.getSize() != nVar) {
			cerr<<"initializing the target"<<endl;
		}
		if (x.getSize() != obj.getSize()){
			cerr<<"X:"<<x.getSize()<<" obj: "<<obj.getSize()<<endl;
		}

		//create objective as 0
		model.add(IloMaximize(env, 0));

		IloRangeArray con(env);
		//add p constraints
		for (int j = 0; j < n; j ++){
			IloExpr pcons(env);
			for(int i = 0; i < m;  i ++){
				pcons += x[j*m + i] * 1;
			}
			con.add(M.P[j] <= pcons <= M.P[j]);
		}
		//cout<<"p constraints added"<<endl;
		//add R constraints, e.g. c00 constraints
		int base = m*n;
		for(int i = 0; i < nPair; i ++){
			int begin = m*i;
			int end = m*(i+1)-1;

			IloExpr rcons(env);
			for(int j = begin; j <= end; j ++){
				rcons += x[base + j] * 1;
			}
			con.add(M.LR[i] <= rcons <= M.LR[i]);
		}
		//cout<<"r constraints added"<<endl;

		//add x_ij ~ axu_ijk constraints
		for (int j = 0; j < n-1; j ++) {
			for (int k = j+1; k < n; k ++ ) {
				for (int i = 0; i < m; i ++ ) {
					int d1 = j*m + i;
					int d2 = k*m + i;
					int dd_relative = calcdd(i, j, k, m, n);
					int dd = dd_relative + m*n - 1;
					//cout<<d1<<" "<<d2<<" "<<dd_relative<<"/"<<dd<<endl;
					IloExpr auxConsTop(env);
					auxConsTop = - x[d1] - x[d2] + 2 * x[dd];
					con.add(-INFINITY <= auxConsTop <= 0);

					IloExpr auxConsBtn(env);
					auxConsBtn =   x[d1] + x[d2] - 2 * x[dd];
					con.add(-INFINITY <= auxConsBtn <= 1);
				}
			}
		}

		//adding fixing constraints for each row
		for(int f = 0; f < m; f++){
			int fixVarBeginIndex = m*n + nAuxVar + m*n*f;	//beginning of aux variable
			int *fixVar = M.M[f];
			for(int i = 0; i < m; i ++){		//for each row
				for(int j = 0; j < n; j ++){	//for each col
					int x_ij = j*m+i;	//index variable
					int fixVarIdx = fixVarBeginIndex + j*m+i;

					IloExpr fixVarCons1(env);
					fixVarCons1 = x[fixVarIdx] - x[x_ij];
					con.add(-INFINITY <= fixVarCons1 <= fixVar[j]);

					IloExpr fixVarCons2(env);
					fixVarCons2 = x[fixVarIdx] + x[x_ij];
					con.add(-INFINITY <= fixVarCons2 <= 2 - fixVar[j]);

					IloExpr fixVarCons3(env);
					fixVarCons3 = x[fixVarIdx] + x[x_ij];
					con.add(fixVar[j] <= fixVarCons3 <= INFINITY);

					IloExpr fixVarCons4(env);
					fixVarCons4 = x[fixVarIdx] - x[x_ij];
					con.add(-fixVar[j] <= fixVarCons4 <= INFINITY);
				}
			}
		}

		//constraints for fixing variable
		for(int f = 0; f < m; f++){
			int fixVarBeginIndex = m*n + nAuxVar + f*m*n;
			for(int i = 0; i < m; i ++){
				IloExpr fcons(env);
				for(int j = 0; j < n; j ++){
					int idx = fixVarBeginIndex + j*m+i;
					fcons += x[idx]*1;
				}
				con.add(1 <= fcons <= n);	//for each row, at least one bit different.
			}
		}

		model.add(con);

		IloCplex cplex(model);
		cplex.exportModel("exclude.lp");

		//set parameters before populate solutions
//		cplex.setParam(IloCplex::SolnPoolIntensity, 4);	//enum all solutions
//		cplex.setParam(IloCplex::SolnPoolAGap, 0.0);
//		cplex.setParam(IloCplex::PopulateLim, 1000);	//the top of solutions
//		cplex.setParam(IloCplex::PopulateSolLim, );
//		cplex.setParam(IloCplex::SolnPoolCapacity, 10000000000);
		bool found = cplex.solve();

		if(found){
			cout<<"found solution"<<endl;
			IloNumArray vals(env);
			cplex.getValues(vals, x);
			sp->addToPool(vals, m, n);
			c2<<"original"<<endl;
			M.printMatrix(c2);
			c2<<"sln"<<endl;
			sp->printPool(c2);
		}else{
			env.out()<<"sln status"<<cplex.getStatus()<<endl;
			cout<<"fail"<<endl;
		}
		//print solutions
//		//check the solutions
//		numSolution = cplex.getSolnPoolNsolns();
//		int nRealSln = 0;
//		//TODO access solutions and remove duplicate solutions
//		for(int si = 0; si < numSolution; si++){
//			IloNumArray vals(env);
//			cplex.getValues(vals, x, si);
//			if(sp->addToPool(vals, m, n)){
//				nRealSln ++;
//			}
//		}

		//sp->printPool(cout);

		env.out()<<"the solution pool contains "<<numSolution<<" solutions"<<endl;
		env.end();
	}
	catch (IloException &e)
	{
		cerr<<"concert exception caught: "<< e << endl;
	}

	return sp;
}

