/*
 * AuxFunc.cpp
 *
 *  Created on: Sep 17, 2010
 *      Author: xzhou
 */

#include "AuxFunc.h"
#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>

using namespace std;

int ** readMatrixFromFile(string fileName, int &nInd, int &nSnp){
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
	int **M = new int*[nInd];
	for(int i = 0; i < nInd; i++) {
		M[i] = new int[nSnp];
	}

	snpFile.close();
	ifstream snpFile1(fileName.c_str());
	//read the file to the array
	for(int i = 0; i < nInd; i ++) {
		for(int j = 0; j < nSnp; j ++ ) {
			snpFile1>>M[i][j];
			//cout<<M[i][j]<<" ";
		}
		//cout<<endl;
	}
	snpFile.close();
	return M;
}

//allocate m by n matrix, initilized it with defaultVal
int ** newInt2d(int m, int n, int defaultVal) {
	int **tmpM = new int*[m];
	for(int i = 0; i < m; i++) {
		tmpM[i] = new int[n];
	}
	//initialize
	for(int i = 0; i < m; i ++) {
		for(int j = 0; j < n; j ++){
			tmpM[i][j] = defaultVal;
		}
	}
	return tmpM;
}

//calculate n choose k
long nchoosek(int n, int k){
	if (k == 2) {
		return n*(n-1)/2;
	}
	else {
		// TODO implement this
		cerr<<"not implemented yet"<<endl;
		assert(false);
		return 0;
	}
}

int calcPairIndex(int i, int j, int n){
	int d = 0;
	if(i > j){
		cerr<<"i < j faile"<<endl;
		assert(false);
	}
	d = i*(2*n-i-1)/2;
	d = d + j - i;
	return d;
}

int calcdd(int i, int j, int k, int m, int n){
	int pairIndex = calcPairIndex(j, k, n);
	return (pairIndex-1)*m + i + 1;
}

//return a matrix of size subm x subn
int ** randSubMatrix(int **M, int m, int n, int subm, int subn){
	//TODO need test case
	int ** retM = newInt2d(subm, subn);

	if(subm > m || subn > n){
		cerr<<"can not get a sub matrix larger than the original one"<<endl;
		return NULL;
	}

	int mStart = rand()%(m-subm);
	int nStart = rand()%(n-subn);

	cout<<mStart<<":"<<nStart<<endl;
	cout<<subm+mStart<<":"<<subn+nStart<<endl;

	for(int i = 0; i < subm; i ++ ){
		for(int j = 0; j < subn; j ++){
			retM[i][j] = M[i+mStart][j+nStart];
		}
	}
	return retM;
}

void delete2d(int **M, int m, int n){
	for(int i = 0; i < m; i ++){
		delete[] M[i];
	}
	delete[] M;
}

void print2d(int **M, int m, int n){
	for(int i = 0; i < m; i ++){
		for(int j = 0; j < n; j ++){
			printf("%1d ", M[i][j]);
		}
		printf("\n");
	}
}
