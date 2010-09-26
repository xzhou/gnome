/*
 * AuxFunc.cpp
 *
 *  Created on: Sep 17, 2010
 *      Author: xzhou
 */

#include "AuxFunc.h"
#include <assert.h>
#include <iostream>

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


