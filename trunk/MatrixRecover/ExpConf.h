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
	static int nMin = 10;	//max snps
	static int nMax = 20;	//min snps
	static int mMin = 10;	//min individuals
	static int mMax = 20;	//max individuals
	ExpConf();
	virtual ~ExpConf();
};

#endif /* EXPCONF_H_ */
