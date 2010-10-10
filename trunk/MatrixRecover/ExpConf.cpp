/*
 * ExpConf.cpp
 *
 *  Created on: Oct 1, 2010
 *      Author: xzhou
 */

#include "ExpConf.h"
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

ExpConf::ExpConf(string fileName) {
	// TODO Auto-generated constructor stub
//	nMin = 9;
//	nMax = 30;
//	mMin = 9;
//	mMax = 10;
//	repeat = 10;
//	diff = 5;
	readConf(fileName);
}

ExpConf::~ExpConf() {
	// TODO Auto-generated destructor stub
}

void ExpConf::readConf(string fileName){
	ifstream f(fileName.c_str());
	string line;
	getline(f, line);//ignore the header
	getline(f, line);
	stringstream s(line);
	s>>mMin>>mMax>>nMin>>nMax>>repeat;
}

string ExpConf::toString(){
	stringstream s;
	s<<mMin<<"-"<<mMax<<":"<<nMin<<"-"<<nMax<<"x"<<repeat;
	return s.str();
}
