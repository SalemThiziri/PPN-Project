/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */ 
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <mpi.h>
#include "blaslapack_headers.h"


struct func_args_s
{   int ntot;
    int rank;
    int nproc;
    int nloc;
    int istart;
    int iend;
    double* xn;
    double* yn;
    double*yloc;
    double *A;
   
};
typedef struct func_args_s func_args_t;



double* AllouerTableau(int taille);
void set_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la);
void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int *la, int *kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int *la, int *kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);

//void produitParallel(double * A,double *x,int istart, int iend,int nloc,int ntot, int lab, int kl, int ku, double *yn,double *xn, double *y,int nproc, int rank,double * yloc);
void produitParallel( double * x,double *Y,int lab, int kl, int ku,func_args_t *all_args);
void  set_dense_X0_1D(double* X0, int* la);
double*  GradientConjuge(double*AB,double*X0 ,double* tableau,int *ku,int *kl,int *lab, double epsilon, int la, int ItMax);

//double* GradientConjugeParallel(double *A,int istart, int iend,int nloc,double *Y, double *yn,double * xn, double*X0 ,double* tableau,int *ku,int *kl,int *lab,  double epsilon, int la,int ItMax,int nproc,int rank,double *yloc);
double* GradientConjugeParallel( double*X0 ,double* tableau,int *ku,int *kl,int *lab,  double epsilon, int la,int ItMax,func_args_t *all_args);
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);
void eig_poisson1D(double* eigval, int *la);
double eigmax_poisson1D(int *la);
double eigmin_poisson1D(int *la);
double richardson_alpha_opt(int *la);
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit);
