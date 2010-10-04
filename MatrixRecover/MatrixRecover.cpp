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

	//the total size of a snp file

	SnpMatrix* allMatrix = readMatrixFromFile(fileName);
	Solver slv;
	ExpConf conf;

	//random sample a matrix of size m and n and solve it, we need to explore the space
	outputFile<<"m\tn\tsn\n"<<endl;
	for(int m= conf.mMin; m < conf.mMax; m ++){
		int nBase = 2*m/log(m+1);
		for(int k = 0; k < 2*conf.diff; k++){
			int n = nBase + k - conf.diff;
			if(n < 0) {
				n = 1;
			}
			for(int r = 0; r < conf.repeat; r++){
				SnpMatrix *subM = allMatrix->randSubMatrix(m, n);
				int sn1 = slv.solve(*subM);
				outputFile<<m<<"\t"<<n<<"\t"<<sn1<<endl;
				delete subM;
			}
		}
	}
	delete allMatrix;
}
