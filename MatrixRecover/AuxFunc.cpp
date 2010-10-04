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

void trimSpace(string& str)
{
    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if(( string::npos == startpos ) || ( string::npos == endpos))
    {
        str = "";
    }
    else
    {
        str = str.substr( startpos, endpos-startpos+1 );
    }
}
