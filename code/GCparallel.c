/******************************************/
/*               */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include<mpi.h>
#include"time.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  MPI_Init(&argc, &argv);

 
  int nbpoints, la;
  int ku, kl, kv, lab;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X,*X0;
  double temp, relres;

  NRHS=1;
  nbpoints=atoi(argv[1]);
  la=nbpoints-2;
  T0=5.0;
  T1=0.0;
  
  func_args_t  *all_args;
  all_args=(func_args_t*)malloc(sizeof(func_args_t));
 
	int nloc, istart, iend;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &all_args->rank);
	MPI_Comm_size(MPI_COMM_WORLD, &all_args->nproc);
if(all_args->rank==0)
  printf("--------- Poisson 1D with GC---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  X0=(double *) malloc(sizeof(double)*la);
  
  
  //la décomposition de domain
  set_grid_points_1D(X, &la);
  //le vecteur second membre
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  //la solution analystique 
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  // la solution initiale 
  set_dense_X0_1D(X0,&la);
  
   
     // ecrire le vecteur x0
     write_vec(X0, &la, "X0.dat");
     //ecrire le le vecteur second membre 
     //write_vec(RHS, &la, "RHS.dat");
     //ecrire le vecteur solution
     write_vec(EX_SOL, &la, "EX_SOL.dat");
     //ecrire le decoposition
     //write_vec(X, &la, "X_grid.dat");
     
  double epsilon=0.000001;
  kv=0;//pour enlever la colonne des 0
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
 //double deb,fin;

 double* x_sol=(double *) malloc(sizeof(double)*la);
 
 
 
 
	
	int q = la/all_args->nproc;
	int r = la%all_args->nproc;
        all_args->ntot=la;
	
	if(all_args->rank<r){
		all_args->nloc = q+1;
		all_args->istart = all_args->rank * (q+1);
	}else{
		all_args->nloc = q;
		all_args->istart = r * (q+1) + (all_args->rank-r) * q;
	}
	all_args->iend = all_args->istart+all_args->nloc;
	
	all_args->xn = AllouerTableau(all_args->nloc);

	
	all_args->A = AllouerTableau(all_args->nloc*lab);	
	set_GB_operator_colMajor_poisson1D(all_args->A, &lab, &all_args->nloc, &kv);
	 //double*     all_args->Y=(double*)calloc(la,sizeof(double));
	
	
	//write_GB_operator_colMajor_poisson1D(A, &lab, &nloc, "AB.dat");
	
	all_args->yn = AllouerTableau(all_args->nloc);
 
        all_args->yloc = (double*)calloc(la, sizeof(double));
 
 
 
  //x_sol=GradientConjugeParallel(A,istart,iend,nloc,Y,yn,xn,X0 ,RHS,&ku,&kl,&lab,epsilon,la, 400,nproc,rank,yloc);
 
 
  x_sol=GradientConjugeParallel(X0 ,RHS,&ku,&kl,&lab,epsilon,la, la, all_args);
  

 if(all_args->rank==0){
 write_vec(x_sol, &la, "x_sol.dat");

    printf("\n\n--------- End -----------\n");

}
  free(RHS);
  free(EX_SOL);
  free(X);
  free(X0);
  free(x_sol);
 free(all_args->xn);
 free(all_args->yn);
 free(all_args->A);
 free(all_args->yloc);


  MPI_Finalize();


}
