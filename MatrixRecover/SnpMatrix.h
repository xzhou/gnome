/*
 * SnpMatrix.h
 *	This is a matrix representation of the snp matrix.
 *  Created on: Sep 12, 2010
 *      Author: xzhou
 */

#ifndef SNPMATRIX_H_
#define SNPMATRIX_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "AuxFunc.h"

using namespace std;

class SnpMatrix {
private:

public:
	int **M;	//binary matrix
	int nInd;	//number of individuals, aka number rows of M
	int nSnp;	//number of SNPs, aka number of cols of M
	int *P;		//single allele frequency
	int **R;	//pairwise, r
	int *LR;	//linear representation of R;


	SnpMatrix(string fileName);		//read matrix from a file
	SnpMatrix(int **M, int m, int n);
	SnpMatrix(int m, int n); //random matrix generator
	virtual ~SnpMatrix();
	void readMatrixFromFile(string fileName);
	void calculateP();	//calculate single allele freq
	void calculateR();	//calculate pairwise allele freq
	//int ** newInt2d(int m, int n);	//create 2d array of m row, n col
	void clearMem();	//clear memories
	void convLinearR();	//convert to one dimension pairwise R
	SnpMatrix* randSubMatrix(int subm, int subn);
	void printMatrix(ostream &o = cout);
	void calculate();
};

SnpMatrix* readMatrixFromFile(string fileName);

#endif /* SNPMATRIX_H_ */
