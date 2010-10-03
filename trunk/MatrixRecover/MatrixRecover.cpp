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
	for(int mi = conf.mMin; mi < conf.mMax; mi ++){
		for(int j = 0; j < conf.repeat; j++){
			int mj_1 =	2*mi/log(mi+1);
			int **M1 = randSubMatrix(allMatrix, mAll, nAll, mi, mj_1);
			int sn1 = slv.solve(M1, mi, mj_1);
			//ajacent scale matrix
			int mj_2 = mj_1 + 1;
			int mj_3 = mj_1 - 1;
			int **M2 = randSubMatrix(allMatrix, mAll, nAll, mi, mj_2);
			int **M3 = randSubMatrix(allMatrix, mAll, nAll, mi, mj_3);
			int sn2 = slv.solve(M2, mi, mj_2);
			int sn3 = slv.solve(M3, mi, mj_3);

			//cout<<"mi = "<<mi<<" mj = "<<mj<<" mj_2 = "<<mj_2<<" mj_3 = "<<mj_3<<endl;
			outputFile<<mi<<" "<< mj_1<< " "<<mj_2<<" "<<mj_3<<endl;
			outputFile<<sn1<<" "<< sn2<<" "<<sn3<<endl<<endl;;
			delete2d(M1, mi, mj_1);
			delete2d(M2, mi, mj_2);
			delete2d(M3, mi, mj_3);
		}
	}

	if(allMatrix) delete2d(allMatrix, mAll, nAll);
}
