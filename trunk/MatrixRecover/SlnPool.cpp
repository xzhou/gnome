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
//#include <openssl/md5.h>
#include <stdio.h>
#include <sstream>

SlnPool::SlnPool() {
	// TODO Auto-generated constructor stub

}

SlnPool::~SlnPool() {
	// TODO Auto-generated destructor stub
}

//try to add the new solution to the pool, return 1 if new solution found
//return 0 if existed
int SlnPool::addToPool(const int *M, int m, int n) {
	string md = m2s(M, m, n);
	return addToPool(md);
}

int SlnPool::addToPool(IloNumArray &M, int m, int n){
	string s = m2s(M, m, n);
	return addToPool(s);
}

int SlnPool::addToPool(string md){
	if(solutionPool[md] == 0){
		solutionPool[md] ++;
		cout<<md<<endl;
		return 1;
	}else{
		solutionPool[md] ++;
		return 0;
	}
}

bool SlnPool::existInPool(string s) {
	SLNMAP_ITR it;
	it = solutionPool.find(s);
	if(it == solutionPool.end()) return false;
	return true;
}

string SlnPool::m2s(IloNumArray &M, int m, int n){
	//convert a matrix to a	string
	//convert
	vector<string> ss(m);
	for(int i = 0; i < m; i ++){
		stringstream s;
		for(int j = 0; j < n; j++)
		{
			s<<M[j*m+i];
//			cout<<M[j*m+i];
		}
		s<<endl;
		ss[i] = s.str();
	}
	sort(ss.begin(), ss.end());
	//return matrix s
	string	ms;
	for(int i = 0; i < m; i ++ ) {
		ms += ss[i];
		cout<<ss[i]<<endl;
	}
	cout << ms << endl;
	return ms;
}

string SlnPool::m2s(const int *M, int m, int n){
	//convert a matrix to a	string
	//convert
	vector<string> ss(m);
	for(int i = 0; i < m; i ++){
		stringstream s;
		for(int j = 0; j < n; j++)
		{
			s<<M[j*m+i];
//			cout<<M[j*m+i];
		}
		s<<endl;
		ss[i] = s.str();
	}
	sort(ss.begin(), ss.end());
	//return matrix s
	string	ms;
	for(int i = 0; i < m; i ++ ) {
		ms += ss[i];
		cout<<ss[i]<<endl;
	}
	cout << ms << endl;
	return ms;
}

string SlnPool::s2signature(string s) {
	//calculate the hash value of a string
	//cout<<s<<endl;

//	unsigned char result[MD5_DIGEST_LENGTH];
//	unsigned char *md =  MD5((unsigned char *)s.c_str(), s.length(), result);
//	char tmp[2];
	string ss = "";
//	for(int i = 0; i < MD5_DIGEST_LENGTH; i++){
//		sprintf(tmp, "%02x", result[i]);
//		printf("%02x", result[i]);
//		ss += tmp;
//	}
	cout<<ss;
	return ss;
}


