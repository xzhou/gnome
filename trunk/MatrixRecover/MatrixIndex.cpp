/*
 * MatrixIndex.cpp
 *
 *  Created on: Sep 26, 2010
 *      Author: xzhou
 */

#include "MatrixIndex.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <openssl/md5.h>

MatrixIndex::MatrixIndex() {
	// TODO Auto-generated constructor stub

}

MatrixIndex::~MatrixIndex() {
	// TODO Auto-generated destructor stub
}

int MatrixIndex::addToPool(int *M, int m, int n) {
	return 0;
}

bool MatrixIndex::existInPool(int *M, int m, int n) {
	return 0;
}

string MatrixIndex::m2s(int *M, int m, int n){
	//convert a matrix to a	string
	//convert
	//sort
	char s[m*n];
	for(int i = 0; i < m*n; i ++){
		s[i] = M[i] + 48;
	}
	string longS = string(s);
	vector<string> ss(m);
	for(int i = 0; i < n; i ++){
		ss[i] = longS.substr(i*m, (i+1)*m-1);
	}
	sort(ss.begin(), ss.end());

	//return matrix s
	string	ms;
	for(int i = 0; i < m; i ++ ) {
		ms += ss[i];
	}
	return ms;
}

string MatrixIndex::s2signature(string s) {
	//calculate the hash value of a string
	unsigned char result[MD5_DIGEST_LENGTH];
	unsigned char *md =  MD5(s.c_str(), s.length(), result);
	return "";
}
