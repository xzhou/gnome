/*
 * MatrixIndex.cpp
 *
 *  Created on: Sep 26, 2010
 *      Author: xzhou
 */

#include "SlnPool.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <openssl/md5.h>

SlnPool::SlnPool() {
	// TODO Auto-generated constructor stub

}

SlnPool::~SlnPool() {
	// TODO Auto-generated destructor stub
}

int SlnPool::addToPool(int *M, int m, int n) {
	string md = m2s(M, m, n);
	bool exist = existInPool(md);
	if (! exist) {
		solutionPool.insert(md);
		return true;
	}
	return false;
}


bool SlnPool::existInPool(string s) {

	set<string>::iterator it;
	it = solutionPool.find(s);
	if(it == solutionPool.end()) return false;
	return true;
}

string SlnPool::m2s(int *M, int m, int n){
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

	return s2signature(ms);
}

string SlnPool::s2signature(string s) {
	//calculate the hash value of a string
	unsigned char result[MD5_DIGEST_LENGTH+1];
	unsigned char *md =  MD5((unsigned char *)s.c_str(), s.length(), result);
	char tmp[2];
	string ss = "";
	for(int i = 0; i < MD5_DIGEST_LENGTH; i++){
		sprintf(tmp, "%02x", result[i]);
		ss += tmp;
	}
	return ss;
}


