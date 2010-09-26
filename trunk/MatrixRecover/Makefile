#============================================================================
# Name        : Makefile
# Author      : Xiaoyong Zhou
# Version     :
# Copyright   : All Rights Reserved @ zhou@indiana.edu
# Description : GNU makefile
#============================================================================

SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /home/xzhou/Applications/ILOG/CPLEX_Studio_AcademicResearch122/cplex
CONCERTDIR    = /home/xzhou/Applications/ILOG/CPLEX_Studio_AcademicResearch122/concert
# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------
#CCOPT = -m64 -O -fPIC -fexceptions  -DIL_STD
CCOPT = -m64 -O2 -g -Wall -fPIC -fexceptions -DNDEBUG -DIL_STD
COPT  = -m64 -fPIC 

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread
JAVA      = java -d64 -Djava.library.path=$(CPLEXDIR)/bin/x86-64_sles10_4.1 -classpath $(CPLEXJARDIR):


all:
	make all_cpp

execute: all
	make execute_cpp
	
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCC        = $(EXDIR)/src/c
EXSRCCPP      = $(EXDIR)/src/cpp
EXSRCJAVA     = $(EXDIR)/src/java

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------

OBJS = MatrixRecover.o SnpMatrix.o AuxFunc.o
TARGET = MatrixRecover

CPP_EX = $(TARGET)
all_cpp: $(CPP_EX)
execute_cpp: $(CPP_EX)
	./MatrixRecover

# --------make objects--------------------------
%.o: %.cpp
	$(CCC) -c $(CCFLAGS) -o $@ $<

# --------make target---------------------------
$(TARGET): $(OBJS)
	$(CCC) $(CCFLAGS) -o $(TARGET) $(OBJS) $(CCLNFLAGS)
	
#-------- make clean ---------------------------
clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp	
# Local Variables:
# mode: makefile
# End: