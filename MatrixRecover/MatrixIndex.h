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

class MatrixIndex {
	list<string> slnPool;	//a MD5 solution list
public:
	MatrixIndex();
	int addToPool(int **M, int m, int n);
	bool existInPool(int **M, int m, int n);
	virtual ~MatrixIndex();
};

#endif /* MATRIXINDEX_H_ */
