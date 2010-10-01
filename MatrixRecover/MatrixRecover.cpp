//============================================================================
// Name        : SnpSolver.cpp
// Author      : Xiaoyong Zhou
// Version     :
// Copyright   : All Rights Reserved @ zhou@indiana.edu
// Description : MatrixRecover in C, Ansi-style, use IBM cplex to reverse an
//				 binary matrix
//============================================================================

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "SnpMatrix.h"
#include "AuxFunc.h"
#include "SolutionFilterCallback.h"

ILOSTLBEGIN

// TODO
//we need a database to remember all unique solutions seen so far

int main(int argc, char *argv[]) {
	//read data from file
	string fileName = "";

	int mLim = 0;
	int nLim = 0;

	if(argc == 2){
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

	outputFile = ofstream(outputFileName);

	//the total size of a snp file
	int **allMatrix = NULL;
	int mAll = 0;
	int nAll = 0;

	readMatrixFromFile(fileName.c_str(), allMatrix, mAll, nAll);

	//random sample a matrix of size m and n and solve it, we need to explore the space

	int **M = randSubMatrix(allMatrix, mAll, nAll, 10, 10);

	ExpConf conf();


	if(allMatrix) delete allMatrix;
}
