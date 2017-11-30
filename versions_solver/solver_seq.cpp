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
  Vector residu      = b - A*u;
  double residuNorm0 = 0;
  for (int i = 0 ; i < m.nbOfNodes ; i++) residuNorm0 += pow(residu(i),2);
  double residuNorm  = tol*residuNorm0 + 1;  // au cas où le tol initial soit négatif ...
  int it = 0;
  while (residuNorm > tol * residuNorm0 && it < maxit){
    
    // Compute N*u
    Vector Nu = N*u;
    exchangeAddInterfMPI(Nu, m);
    
    // Update field
    for(int n=0; n<m.nbOfNodes; n++)
      u(n) = 1/Mdiag(n) * (Nu(n) + b(n));

    
    // Update residual and iterator
    residu = b - A*u;
    residuNorm = 0;
    for (int i=0; i<m.nbOfNodes; i++) residuNorm += pow(residu(i),2);
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
}

