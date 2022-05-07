

 ____ __ ___ ____ ________ _________ ______   __ ________ ___________ ______ __ _______
|                                                                                      |
|*/Parallel implementation of an iterative algorithm for conjugate gradient method./*  |                                                                               |
|____ __ ___ ____ ________ _________ ______   __ ________ ___________ ______ __ _______|





-----------------------------------/****Parallel Computing Project*****/--------------------------------



----------------------------*****Notes about our implementation*****--------------------------------------

  -->  We’ve created one library  lib_poisson1D.c (which contains main implementation) 
    
  -->  We’ve created 2 driver programs :
  
     * GC.c : the sequential  version  solves the given  problem .
        
     * GCparallel.c: the parallel version .
       
 ---> The inputs are: The matrix size.
 
 ---> The  outputs are :
      
  * AB.dat:  represent  the storage of the matrix with col-major.
   
  * X_grid.dat: the  discretizationof the interval of the resoulotion.
  
  * RHS.dat:The second membre vector .
  
  * X0.dat: The initial solution .
  
  * EX_SOL.dat: The analytical solution.
  
  * x_sol.dat: The approximat solution calculated by the GC program .
  
  
--------------------------------------*****Compilation*****-----------------------------------------------------
        
       To compile the code use the following command:

  $ make

To clean binaries:

  $ make clean


--------------------------------------***** Execution*****-----------------------------------------------------

	To execute the sequential  version:      ./bin/GC size_of_matrix
	To execute the parallel version  :       ./bin/GCparallel size_of_matrix
	
