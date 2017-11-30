#include "headers.hpp"

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
  Vector Mdiag(A.rows());
  SpMatrix N(A.rows(), A.cols());
  for(int k = 0; k < A.outerSize(); ++k){
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, m);
  
  // Jacobi solver
  double residuNorm = tol*2;
  int it = 0;
  while (residuNorm > tol && it < maxit){
    
    // Compute N*u
    Vector Nu = N*u;
    exchangeAddInterfMPI(Nu, m);
    
    // Update field
    for(int n=0; n<m.nbOfNodes; n++)
      u(n) = 1/Mdiag(n) * (Nu(n) + b(n));
    
    // Update residual and iterator
    if((it % 10) == 0){
      if(myRank == 0)
        printf("\r   %i %e", it, residuNorm);
    }
    it++;
  }
  
  if(myRank == 0){
    printf("\r   -> final iteration: %i (prescribed max: %i)\n", it, maxit);
    printf("   -> final residual: %e (prescribed tol: %e)\n", residuNorm, tol);
  }
}

