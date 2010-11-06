/*
 * Experiments.cpp
 *
 *  Created on: Nov 5, 2010
 *      Author: xzhou
 */

#include "Experiments.h"
#include <ilcplex/ilocplex.h>
#include "ExpConf.h"
#include <fstream>


Experiments::Experiments() {
	// TODO Auto-generated constructor stub

	conf = ExpConf("expconf");
	cout<<conf.toString()<<endl;
	string confstr = conf.toString();

	//read data from file
	string snpFileName = "test80.txt";

	string outputFileName = "slnct_"+confstr+".log";
	slnCtLog = new ofstream(outputFileName.c_str());

	string name = "sln_v_"+confstr+".log";
	slnLog = new ofstream(name.c_str());//verbose log

}


Experiments::~Experiments() {

}

void Experiments::findMultipleSolution(){
	Solver slv;
	*slnCtLog<<"m\tn\tsn\n"<<endl;
	for(int m = conf.mMin; m <= conf.mMax; m ++){
		int nBase = 2*m/log(m+1);
		int N[3];
		N[0] = nBase/2;	//smaller space
		N[1] = nBase;	//border space
		N[2] = nBase*2;	//"unique" space
		for(int n = conf.nMin; n <= conf.nMax; n++){
			//int n = nBase + k - conf.diff;
			//int n = N[k];

			if(n < 0) {
				n = 1;
			}
			for(int r = 0; r < conf.repeat; r++){
				//SnpMatrix *subM = allMatrix->randSubMatrix(m, n);		//use human genome to generate matrix
				//SlnPool* sp1 = slv.solve(*subM);
				SnpMatrix *subM = new SnpMatrix(m, n);	//use
				SlnPool* sp2 = slv.solveAndFilter(*subM);
				stringstream ss;
				ss<<m<<"\t"<<n<<"\t"<<sp2->getNumSln()<<endl;
				string tmp = ss.str();
				cout<<tmp;
				*slnCtLog<<tmp<<endl;
				*slnLog<<"**********"<<endl;
				*slnLog<<tmp;
				subM->printMatrix(*slnLog);
				*slnLog<<"--"<<endl;
				sp2->printPool(*slnLog);
				delete subM;
			}
		}
	}
}

void Experiments::largeScaleFixRow(){
	//this method generate large matrix or sample large matrix from real genome sequence
	Solver slv;
	int m = conf.mMin;
	int n = conf.nMin;
	int repeat = conf.repeat;
	for(int i = 0; i < repeat; i ++){
		SnpMatrix *M = new SnpMatrix(m, n);
		for(int j = 0; j < m; j ++){
			SlnPool* sp = slv.fixAndSolve(*M, j);
			int nSln = sp->getNumSln();
			*slnCtLog<<"Sln = "<<nSln<<endl;
			if(nSln == 0){
				*slnCtLog<<"fixed: j = "<<j<<endl;
			}
		}
	}
}

void samllScaleFixRow(int m, int n, int repeat){

}

