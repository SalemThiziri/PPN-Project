
############################################
#########    makefile        ###############



################################
# Variables for this makefile
################################
# 
CC=gcc
CC2=mpicc

# 
# -- Compiler Option
#
#OPTC=-O3 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
#OPTC=-O2 -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
#OPTC=-O -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
OPTC=-Ofast  -fomit-frame-pointer -fPIC -mavx -DAdd_ -DF77_INTEGER=int -DStringSunStyle
# -- Directories
DIR=.
DIRSRC=$(DIR)/src


# -- librairies
LIBS=-lblas -lm

# -- Include directories
INCLBLAS= -I /usr/include

INCL= -I $(DIR)/include $(INCLBLAS) 


OBJSEQ= lib_poisson1D.o GC.o
OBJPARALLEL= lib_poisson1D.o GCparallel.o
#

all: bin/GC bin/GCparallel

GC: bin/GC
 
GCparallel: bin/GCparallel

lib_poisson1D.o: $(DIRSRC)/lib_poisson1D.c
	$(CC2) $(OPTC) -c $(INCL) $(DIRSRC)/lib_poisson1D.c 

  

GC.o: $(DIRSRC)/GC.c
	$(CC2) $(OPTC) -c $(INCL) $(DIRSRC)/GC.c  
	
GCparallel.o: $(DIRSRC)/GCparallel.c
	$(CC2) $(OPTC) -c $(INCL) $(DIRSRC)/GCparallel.c 

bin/GC: $(OBJSEQ)
	$(CC2) -o bin/GC $(OPTC) $(OBJSEQ) $(LIBS)
	
bin/GCparallel: $(OBJPARALLEL)
	$(CC2) -o bin/GCparallel $(OPTC) $(OBJPARALLEL) $(LIBS)	

run_GC:
	bin/GC
	
run_GCparallel:
	bin/GCparallel
	

clean:
	rm *.o bin/*
