/*
 * ExpConf.h
 *
 *  Created on: Oct 1, 2010
 *      Author: xzhou The configuration class for scheduling the experiments
 */

#ifndef EXPCONF_H_
#define EXPCONF_H_

class ExpConf {

public:
	int nMin;	//max snps
	int nMax;	//min snps
	int mMin;	//min individuals
	int mMax;	//max individuals
	ExpConf();
	virtual ~ExpConf();
};

#endif /* EXPCONF_H_ */
