/*
 * MatrixIndex.h
 *
 *  Created on: Sep 26, 2010
 *      Author: xzhou
 *      This class is used to convert an matrix to it's index, the index is a set of string and it's count
 *
 */

#ifndef MATRIXINDEX_H_
#define MATRIXINDEX_H_

#include <list>
#include <map>
#include <set>
#include <string>
using namespace std;

class MatrixIndex {
	set<string> solutionPool;	//this is a solution MD5 code
public:
	MatrixIndex();
	int addToPool(int *M, int m, int n);
	bool existInPool(int *M, int m, int n);
	string m2s(int *M, int m, int n);
	string s2signature(string s);
	virtual ~MatrixIndex();
};

#endif /* MATRIXINDEX_H_ */
