//============================================================================
// Name        : SnpSolver.cpp
// Author      : Xiaoyong Zhou
// Version     :
// Copyright   : All Rights Reserved @ zhou@indiana.edu
// Description : MatrixRecover in C, Ansi-style, use IBM cplex to reverse an
//				 binary matrix
//============================================================================

#include <iostream>
#include <fstream>
#include <ilcplex/ilocplex.h>
#include <math.h>
#include "SnpMatrix.h"
#include "AuxFunc.h"
#include "SolutionFilterCallback.h"
#include "ExpConf.h"
#include "Solver.h"
#include "AuxFunc.h"
#include <sstream>

int my_debug = 0;

ILOSTLBEGIN

// TODO
//we need a database to remember all unique solutions seen so far

int main(int argc, char *argv[]) {

//	//test code
//	string x = "   1 1 1       ";
//	trimSpace(x);
//	cout<<x<<" "<<(x.length()+1)/2<<endl;
//	return 0;


	//read data from file
	string fileName = "test80.txt";
	string outputFileName = "solutionCount.log";

	int mLim = 0;
	int nLim = 0;
	if (argc == 1)
	{
		//do nothing
	}
	else if(argc == 2){
		fileName = string(argv[1]);
		cout<<fileName<<endl;
	}else if(argc == 3){
		fileName = string(argv[1]);
		outputFileName = string(argv[2]);	//the write the result to a file
	}
	else{
		cout<<"usage: ./MatrixRecover [matrix_file_name]"<<endl;
		return 0;
	}

	ofstream outputFile(outputFileName.c_str());
	ofstream slnLog("sln.log");

	Solver slv;
	//test code here
//	SnpMatrix* testMatrix = readMatrixFromFile("test.txt");
//	int n1 = slv.solve(*testMatrix);
//	int n2 = slv.solveAndFilter(*testMatrix);
//	cout<<n1<< "?=" <<n2<<endl;
//	return 0;

	SnpMatrix* allMatrix = readMatrixFromFile(fileName);

	ExpConf conf;
	//random sample a matrix of size m and n and solve it, we need to explore the space
	outputFile<<"m\tn\tsn\n"<<endl;
	for(int m = conf.mMin; m <= conf.mMax; m ++){
		int nBase = 2*m/log(m+1);
		for(int k = 0; k <= 2*conf.diff; k++){
			int n = nBase + k - conf.diff;
			if(n < 0) {
				n = 1;
			}
			for(int r = 0; r < conf.repeat; r++){
				SnpMatrix *subM = allMatrix->randSubMatrix(m, n);
				//SlnPool* sp1 = slv.solve(*subM);
				SlnPool* sp2 = slv.solveAndFilter(*subM);
				stringstream ss;
				ss<<m<<"\t"<<n<<"\t"<<sp2->getNumSln()<<endl;
				string tmp = ss.str();
				cout<<tmp;
				outputFile<<tmp<<endl;

				slnLog<<"**********"<<endl;
				slnLog<<tmp;
				subM->printMatrix(slnLog);
				slnLog<<"--"<<endl;
				sp2->printPool(slnLog);
				delete subM;
			}
		}
	}
	delete allMatrix;
	outputFile.close();
	slnLog.close();
}
