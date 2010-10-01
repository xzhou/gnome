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

typedef map<string, int> SLNMAP;
typedef map<string, int>::iterator SLNMAP_ITR;

class SlnPool {
	SLNMAP solutionPool;	//this is a solution MD5 code
public:
	SlnPool();
	int addToPool(const int *M, int m, int n);
	bool existInPool(string s);
	string m2s(const int *M, int m, int n);
	string s2signature(string s);
	
	virtual ~SlnPool();
};

#endif /* MATRIXINDEX_H_ */
