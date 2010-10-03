/*
 * AuxFunc.h
 *
 *  Created on: Sep 17, 2010
 *      Author: xzhou
 */

#ifndef AUXFUNC_H_
#define AUXFUNC_H_

#include <string>

using namespace std;

int ** newInt2d(int m, int n, int defaultVal=0);	//allocate 2d array
long nchoosek(int n, int k);	//n choose k
int calcdd(int i, int j, int k, int m, int n);	//calculate the index
int calcPairIndex(int i, int j, int n);	//	calculate the pair index for p_ij
void readMatrixFromFile(string fileName, int **M, int &m, int &n);
int ** randSubMatrix(int **M, int m, int n, int subm, int subn);

#endif /* AUXFUNC_H_ */
