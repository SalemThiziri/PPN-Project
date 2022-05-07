/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// ecrire la matrice de poisson 1D

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
        fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}


// initialiser la matrice de poisson1D

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}


// le second membre 
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
} 







//vecteur initial

void set_dense_X0_1D(double* RHS,int* la){
  int jj;
  
  for (jj=0;jj<(*la);jj++){
    RHS[jj]=0.0;
  }
}  



// calcul de la solution anlytique 

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  



// discritisation de l'intervalle 
void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}



//écrire un vecteur 
void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  







double* AllouerTableau(int taille)
{
 double *tableau= aligned_alloc(64,taille*sizeof(double));
 //double *tableau= malloc(taille*sizeof(double)); // La fonction malloc inscrit dans notre pointeur l'adresse qui a été reservée.
	
 return tableau;
}



double norm2(double* tableau, int taille)
{
	double c=0;
  for (int i = 0; i < taille ; ++i)
  {
  	c=c+tableau[i]*tableau[i];
  }
  return sqrt(c);

}




double  MulTabTab(double* tableau1, double* tableau2, int taille)
{
  double c=0;

  for (int i = 0; i < taille ; ++i)
  {
  	c=c+tableau1[i]*tableau2[i];
  }
  return c;
}



void   MulaXTab(double* tableau1, double alpha, int taille,double * aXtableau)
{
	
	for (int i = 0; i < taille; ++i)
	{
		aXtableau[i]=alpha*tableau1[i];
	}
	
}
double *Sub_uv(double * tableau1, double * tableau2, int taille)
	{
		double* tableau3=aligned_alloc(64,taille*sizeof(double));
		//tableau3=AllouerTableau(taille);

		for (int i = 0; i < taille; ++i)
		{
			tableau3[i]=tableau1[i]-tableau2[i];
		}
		return tableau3;
	}
	

void  Add_uv(double * tableau1, double * tableau2, int taille,double *tableau3)
	{
		

		for (int i = 0; i < taille; ++i)
		{
			tableau3[i]=tableau1[i]+tableau2[i];
		}
		
	}


	
	
	
	



void produitParallel(double * x,double *Y, int lab, int kl, int ku,func_args_t *all_args)
{
	for(int i=0; i<all_args->nloc; i++) all_args->xn[i] = x[i+all_args->istart];
	cblas_dgbmv(CblasColMajor, CblasNoTrans, all_args->nloc, all_args->nloc, kl, ku, 1.0, all_args->A, lab, all_args->xn, 1, 0.0, all_args->yn, 1);
		
	if(all_args->nproc>1){
		if(all_args->rank==0){
			all_args->yn[all_args->nloc-1] -= x[all_args->iend];
		}else if(all_args->rank==all_args->nproc-1){
			all_args->yn[0] -= x[all_args->istart-1];
		}else{
			all_args->yn[all_args->nloc-1] -= x[all_args->iend];
			all_args->yn[0] -= x[all_args->istart-1];
		}
	}
	
	for(int i=0; i<all_args->nloc; i++) all_args->yloc[i+all_args->istart] = all_args->yn[i];
	
	MPI_Allreduce(all_args->yloc, Y, all_args->ntot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
}




//the main function  of the Gradiant 


double*  GradientConjuge(double*AB,double*X0 ,double* tableau,int *ku,int *kl,int *lab,  double epsilon, int la,int ItMax)
{       double norm;
	int step=0;
	double alpha;
	double beta;
	double err;
	double *tmp=AllouerTableau(la);
	double *tmp2=AllouerTableau(la);
	double *tmp3=AllouerTableau(la);
	double *x=AllouerTableau(la);
	double *r_old=AllouerTableau(la);
	double *r_new=AllouerTableau(la);
	double* p=AllouerTableau(la);
	double *Ap=AllouerTableau(la);
	double *alphaP=AllouerTableau(la);
      double * Y= (double*)calloc(la, sizeof(double));

  
  
       cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,*kl,*ku,1.,AB,*lab,X0,1,0.,Y,1);
        //write_vec(Y, &la, "Y0.dat");

	 r_old=Sub_uv(tableau,Y, la);

	p     =r_old;
	
	norm   =MulTabTab(r_old,r_old,la);
	
	norm=sqrt(norm);
	
	err=norm;
	
	 //printf("erreur1=%e\n",err);
	
	 //write_vec(x, &la, "x.dat");
    
	while(err>epsilon && step<ItMax)

	{         step=step+1;
		  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,*kl,*ku,1.,AB,*lab,p,1,0.,Ap,1);
                  //write_vec(Ap, &la, "Y1.dat");
		   MulaXTab(p,alpha,la,alphaP);
		   alpha  =MulTabTab(r_old,  r_old,  la)   /   MulTabTab(p,Ap, la);
		   MulaXTab(p, alpha ,la,tmp);
		   Add_uv(x    ,    tmp    ,     la,x);
	           //write_vec(x, &la, "x1.dat");
	           MulaXTab(Ap, alpha, la,tmp2);
  		   r_new =Sub_uv(r_old    , tmp2,la);	 
  		   err    =MulTabTab(r_new,r_new,la);//........calculer la norme 2 pour l erreur;
	           err    =sqrt(err);//........calculer la norme 2 pour l erreur;
	          // printf("errreur=%e\n",err);
  		if (err<epsilon)
  		{
  
  			printf("errrur3=%e\n",err);
  			printf("nombre d'iterations=%d \n",step);
  		        free(r_old);
  		        free(r_new); 
	               free(p); 
	               free(Ap); 
	               free(Y); 
	               free(alphaP);
  			return x;
  		}
  		beta= MulTabTab(r_new,r_new,la) /  MulTabTab(r_old, r_old, la);
  		MulaXTab(p, beta, la,tmp3);
  		Add_uv(r_new,tmp3, la,p);
  		r_old=r_new;

  		
			
	}
	printf("Non-Convergence, le nombre d'iterations max est atteint: ");
	               free(r_new); 
	               free(p); 
	               free(Ap); 
	               free(Y); 
	               free(alphaP);
	               free(tmp);
	               free(tmp2);
	               free(tmp3);
return x;	
}




double* GradientConjugeParallel( double*X0 ,double* tableau,int *ku,int *kl,int *lab,  double epsilon, int la,int ItMax,func_args_t *all_args)
{       double norm;
	int step=0;
	double alpha;
	double beta;
	double err;
	double *tmp=AllouerTableau(la);
	double *tmp2=AllouerTableau(la);
	double *tmp3=AllouerTableau(la);
	double *x=AllouerTableau(la);
	double *r_old=AllouerTableau(la);
	double *r_new=AllouerTableau(la);
	double *p=AllouerTableau(la);
	double *Ap=AllouerTableau(la);
	double *alphaP=AllouerTableau(la);
       double * Y=(double*)calloc(la,sizeof(double));
  

        produitParallel(X0,Y,*lab,*kl,*ku,all_args);
        //write_vec(Y, &la, "Y0.dat");

	r_old=Sub_uv(tableau,Y, la);
	p     =r_old;
	norm   =MulTabTab(r_old,r_old,la);
	norm=sqrt(norm);
	err=norm;
       // printf("erreur=%e\n",err);
	
	// write_vec(x, &la, "x.dat");
    
	while(err>epsilon && step<ItMax)
	{
                 step=step+1;
                 
		 produitParallel(p,Ap,*lab,*kl,*ku,all_args);
                // write_vec(Ap, &la, "Y1.dat");
		 MulaXTab(p,alpha,la,alphaP );
		 alpha  =MulTabTab(r_old,  r_old,  la)   /   MulTabTab(p,Ap, la);
		 MulaXTab(p, alpha ,la,tmp);  
		 Add_uv(x    ,  tmp   ,     la,x);
	         write_vec(x, &la, "x1.dat");
	         MulaXTab(Ap, alpha, la,tmp2); 
                 r_new=Sub_uv(r_old    ,   tmp2   ,     la );	 
  		  err    =MulTabTab(r_new,r_new,la);//........calculer la norme 2 pour l erreur;
	         err    =sqrt(err);//........calculer la norme 2 pour l erreur;
	         
  		if (err<epsilon)
  		{
  		
  		
  		if(all_args->rank==0){
  			printf("l'errreur=%e\n",err);
  			printf("nombre d'iterations=%d\n ",step);}
  		       free(r_old); 
	               free(r_new); 
	               free(p); 
	               free(Ap); 
	               free(Y); 
	               free(alphaP); 
	               free(tmp);
	               free(tmp2);
	               free(tmp3);
  			return x;
  		}
  		beta= MulTabTab(r_new,r_new,la) /  MulTabTab(r_old, r_old, la);
  		MulaXTab(p, beta, la,tmp3);
  		Add_uv(r_new, tmp3, la,p);
  		
  		r_old=r_new;

  	
			
	}
	
	  
	//Message de non convergeance.
	printf("Non-Convergence, le nombre d'iterations max est atteint: ");
                       free(r_new); 
	               free(p); 
	               free(Ap); 
	               free(Y); 
	               free(alphaP); 
	               free(tmp);
	               free(tmp2);
	               free(tmp3);
return x;	
}






