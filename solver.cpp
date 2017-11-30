#include "headers.hpp"
#include <unistd.h>

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, Vector& b, Vector& u, Mesh& m, double tol, int maxit)
{
  if(myRank == 0)
    printf("== jacobi\n");
  
  // Compute the solver matrices
  Vector Mdiag(A.rows()); // de taille le nombre de lignes de A
  SpMatrix N(A.rows(), A.cols());
  for(int k = 0; k < A.outerSize(); ++k){  // max(nb_row,nb_col) -> parcourt soit lignes, soit colonnes
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, m);
  
  // Jacobi solver
  Vector Au = A*u;
  exchangeAddInterfMPI(Au, m);
  Vector residu = b - Au;
  //double residuNorm = sqrt(residu.dot(residu.transpose()));
  double residuNorm = ParScalarProd(residu, residu, m);
  
  Vector residu0 = residu;
  //exchangeAddInterfMPI(residu0, m);
  //double residu0Norm = sqrt(residu0.dot(residu0.transpose()));
  double residu0Norm = ParScalarProd(residu0, residu0, m);
  //double residuNorm0 = 1;
  //double residuNorm  = tol*residuNorm0 + 2;  // au cas où le tol initial soit négatif ...
  int it = 0;
  while (residuNorm > tol * residu0Norm && it < maxit){
    
    // Compute N*u
    Vector Nu = N*u;
    exchangeAddInterfMPI(Nu, m);
    
    // Update field
    for(int n=0; n<m.nbOfNodes; n++)
      u(n) = 1/Mdiag(n) * (Nu(n) + b(n));

    
    // Update residual and iterator
    Vector Au = A*u;
    exchangeAddInterfMPI(Au, m);
    residu = b - Au;
    
    //residu = b - A*u;
    //exchangeAddInterfMPI(residu, m);
    //residuNorm = sqrt(residu.dot(residu.transpose()));
    residuNorm = ParScalarProd(residu, residu, m);
    
    if((it % 10) == 0){
      if(myRank == 0)
        printf("\r   %i %e \n", it, residuNorm);
    }
    it++;
  }
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", residuNorm, tol);
  }

  double resnorm_glob;
  MPI_Reduce(&residuNorm, &resnorm_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(myRank == 0) 
    printf("\r Norme residu glob = %e \n", resnorm_glob);
    printf("\r ------------------------------------ \n");
  
}


void ConjGrad(SpMatrix& A, Vector& b, Vector& u, Mesh& m, double tol, int maxit)
{
  if(myRank == 0)
    printf("== conjugate gradient\n");

  // declare variables used in the loop
  double r0Norm2, rNorm2,r_newNorm2, alpha, beta;
  Vector Ap, r_new, p_new, x_new;
  
  // Initialize the CG algorithm at iteration 0
  Vector x = u;
  //exchangeAddInterfMPI(x, m);
  Vector Ax = A*x;
  exchangeAddInterfMPI(Ax, m);
  Vector r  = b - Ax; // b est defini sur les points interieurs et on rajoute les contributions de Ax aux interfaces
  
  Vector p = r;
  Vector r0 = r;
 
  //r0Norm2 = r0.dot(r0.transpose());
  r0Norm2 = ParScalarProd(r0, r0, m); 
  int it = 0;

  r_newNorm2 = r0Norm2;
  do
    {
      //rNorm2 = r.dot(r.transpose());
       
      rNorm2 = ParScalarProd(r, r, m);
      //if(myRank ==0)
      //printf("\r   -> current residu norm: %e\n", sqrt(rNorm2));
      Ap = A*p;
      exchangeAddInterfMPI(Ap, m);

      //alpha = rNorm2 / Ap.dot(p);
      alpha = rNorm2 / ParScalarProd(Ap, p, m);
      //if(myRank == 0)
      //printf("\r   -> alpha: %e\n", alpha);
      
      x_new = x + alpha*p;
      
      r_new = r - alpha*Ap;

      //r_newNorm2 = r_new.dot(r_new.transpose());
      r_newNorm2 = ParScalarProd(r_new, r_new, m);
      //if(myRank == 0)
      //printf("\r   -> new current residu norm: %e\n", sqrt(r_newNorm2));
      beta = r_newNorm2/rNorm2;

      p_new = r_new + beta*p;
      
      x = x_new; r = r_new; p = p_new;
      it ++;
      
    }
  while (sqrt(r_newNorm2) > (tol * sqrt(r0Norm2)) && it < maxit);

  //exchangeAddInterfMPI(x_new, m);
  u = x_new;
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", sqrt(r_newNorm2), tol);
  }
  
  //double rnorm_glob2;
  //MPI_Reduce(&r_newNorm2, &rnorm_glob2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //if(myRank == 0){
  //printf("\r Norme residu glob = %e \n", sqrt(rnorm_glob2));
  // printf("\r ------------------------------------ \n");
  //}
  
}

