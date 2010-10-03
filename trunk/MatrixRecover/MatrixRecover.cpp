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

ILOSTLBEGIN

// TODO
//we need a database to remember all unique solutions seen so far

int main(int argc, char *argv[]) {
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

	//the total size of a snp file
	int **allMatrix = NULL;
	int mAll = 0;
	int nAll = 0;

	allMatrix = readMatrixFromFile(fileName.c_str(), mAll, nAll);
	print2d(allMatrix, mAll, nAll);
	Solver slv;
	ExpConf conf;

	//random sample a matrix of size m and n and solve it, we need to explore the space
	for(int m= conf.mMin; m < conf.mMax; m ++){
		int nBase = 2*m/log(m+1);
		for(int k = 0; k < 2*conf.diff; k++){
			int n = nBase + k - conf.diff;
			if(n < 0) {
				n = 1;
			}
			for(int r = 0; r < conf.repeat; r++){
				int **M1 = randSubMatrix(allMatrix, mAll, nAll, m, n);
				SnpMatrix M(M1, m, n);
				int sn1 = slv.solve(M);
				outputFile<<"m="<<m<<"\tn="<<n<<"\tsn="<<sn1<<endl;
			}
		}
	}

	if(allMatrix) delete2d(allMatrix, mAll, nAll);
}
