#!/usr/bin/python

import sys

ifname = sys.argv[1];
ofname = sys.argv[2];

if ifname == ofname:
    print "warning, ifname == ofname"
    exit();

inputfile = open(ifname, 'r')
outputfile = open(ofname, 'w')

lines = inputfile.readlines()

outputfile.write("m\tn\tsln\n");

for line in lines:
    if len(line) > 1:
        outputfile.write(line)

inputfile.close();
outputfile.close();
