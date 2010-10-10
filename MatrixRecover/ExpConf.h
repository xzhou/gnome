/*
 * ExpConf.h
 *
 *  Created on: Oct 1, 2010
 *      Author: xzhou The configuration class for scheduling the experiments
 */

#ifndef EXPCONF_H_
#define EXPCONF_H_

#include <string>

using namespace std;


class ExpConf {

public:
	int nMin;	//max snps
	int nMax;	//min snps
	int mMin;	//min individuals
	int mMax;	//max individuals
	int repeat;
	int diff;	//difference to f(m);
	ExpConf(string fileName);
	virtual ~ExpConf();
	string toString();
	void readConf(string fileName);
};

#endif /* EXPCONF_H_ */
