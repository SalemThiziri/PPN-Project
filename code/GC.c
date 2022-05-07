/******************************************/
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include"time.h"
int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
 
  int nbpoints, la;
  int ku, kl, kv, lab;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X,*X0;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=atoi(argv[1]);
  la=nbpoints-2;
  T0=5.0;
  T1=0.0;

  printf("--------- Poisson 1D with GC---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  X0=(double *) malloc(sizeof(double)*la);
  
  //decomposition de domain 
  set_grid_points_1D(X, &la);
  //vecteur second membre 
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  // solution initiale  x0
  set_dense_X0_1D(X0,&la);
  // la solution anlytique 
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  
   // ecrire  le vecteur x0
    write_vec(X0, &la, "X0.dat");
  // ecrire le vecteur second membre 
  write_vec(RHS, &la, "RHS.dat");
  //ecrire la solution analytique 
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  //ecrire la decomposition de domain 
  write_vec(X, &la, "X_grid.dat");
  
  double epsilon=0.000001;
  kv=0;//pour enlever la colonne des 0
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

     AB = (double *) malloc(sizeof(double)*lab*la);
     // initialiser la matrice de poisson( stockage bande cilmajor)
     set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
     // ecrire la matrice de poisson 
     write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");
    
    double* x_sol=(double *) malloc(sizeof(double)*la);
   
    x_sol= GradientConjuge(AB,X0 ,RHS,&ku,&kl,&lab,epsilon,la, la);
  
 
    //ecrire la solution calculée par Gradient Conjugué
    write_vec(x_sol, &la, "x_sol.dat");
    

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(X0);
  free(x_sol);

  printf("\n\n--------- End -----------\n");
}
