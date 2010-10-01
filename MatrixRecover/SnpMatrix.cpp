/*
 * SnpMatrix.cpp
 *
 *  Created on: Sep 12, 2010
 *      Author: xzhou
 */

#include "SnpMatrix.h"
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include "AuxFunc.h"

SnpMatrix::SnpMatrix(string fileName) {
	// TODO Auto-generated constructor stub
	M = NULL;
	nInd = 0;
	nSnp = 0;
	P = NULL;
	R = NULL;
	LR = NULL;

	readMatrixFromFile(fileName);
}

//initialized by
SnpMatrix::SnpMatrix(const int **inM, int m , int n){

	M = newInt2d(m, n);

	//we need to copy the data
	for(int i = 0; i < m; i ++){
		for(int j = 0; j < n; j ++){
			M[i][j] = inM[i][j];
		}
	}

	nInd = m;
	nSnp = n;

	calculateP();
	calculateR();
	convLinearR();
}

void SnpMatrix::clearMem(){
	if(M) delete M;
	if(P) delete P;
	if(R) delete R;
	if(LR) delete LR;
}

SnpMatrix::~SnpMatrix() {
	// TODO Auto-generated destructor stub
	clearMem();
}

void SnpMatrix::readMatrixFromFile(string fileName)
{
	//clear memory in case user call thsi function again
	clearMem();

	ifstream snpFile(fileName.c_str());
	string line;
	if(snpFile.is_open()){
		//count the number of lines
		while(getline(snpFile, line)) {
			nInd ++;
			if(nInd == 1){
				nSnp = (line.length()+1)/2;
			}
		}
	}
	cout<<nInd<<"x"<<nSnp<<" matrix read"<<endl;

	//dynamically allocate array;
	M = new int*[nInd];
	for(int i = 0; i < nInd; i++) {
		M[i] = new int[nSnp];
	}

	snpFile.close();
	ifstream snpFile1(fileName.c_str());
	//read the file to the array
	for(int i = 0; i < nInd; i ++) {
		for(int j = 0; j < nSnp; j ++ ) {
			snpFile1>>M[i][j];
			cout<<M[i][j]<<" ";
		}
		cout<<endl;
	}
	snpFile.close();

	//initialize P and R;
	calculateP();
	calculateR();
	convLinearR();
}

void SnpMatrix::calculateP(){
	if(M == NULL){
		cout<<"class not initialized"<<endl;
		return;
	}

	P = new int[nSnp];

	for(int j = 0; j < nSnp; j++){
		int sum = 0;
		for(int i = 0; i < nInd; i ++){
			sum += M[i][j];
		}
		P[j] = sum;
	}
}

void SnpMatrix::calculateR(){
	if(M == NULL){
		cout<<"class not initialized"<<endl;
		return;
	}

	R = newInt2d(nSnp, nSnp);

	for(int i = 0; i < nSnp-1; i++){
		for(int j = i+1; j < nSnp; j++){
			R[i][j] = 0;
			R[j][i] = 0;
			for(int k = 0; k < nInd; k ++){
				if(M[k][i] == 1 && M[k][j] == 1){
					R[i][j] ++;
					R[j][i] ++;
				}
			}
		}
	}
}

//convert it to linear R
void SnpMatrix::convLinearR(){
	int n = nSnp;
	int nPair = nchoosek(n, 2);
	LR = new int[nPair];
	int k = 0;
	for (int i = 0; i < n-1; i ++) {
		for (int j = i+1; j < n; j ++) {
			LR[k++] = R[i][j];
		}
	}
}




