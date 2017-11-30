#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  
  // 1. Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
  
  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  Mesh m;
  readMsh(m, "benchmark/carre_h0.05_N4.msh");
  buildListsNodesMPI(m);
  buildLocalNumbering(m);
  
  // 3. Build problem (fields and system)
  Vector uNum(m.nbOfNodes);
  Vector uExa(m.nbOfNodes);
  Vector f(m.nbOfNodes);
  for(int i=0; i<m.nbOfNodes; ++i){
    double x = m.coords(i,0);
    double y = m.coords(i,1);
    uNum(i) = 0.;
    //uExa(i) = cos(M_PI*x)*cos(2*M_PI*y); // Neumann
    uExa(i) = sin(M_PI*x)*sin(2*M_PI*y); // dirichlet
    f(i) = (1+5*pow(M_PI,2)) * uExa(i);
  }
  
  Problem p;
  double alpha = 1;
  buildLinearSystem(p,m,alpha,f);
  buildDirichletBC(p,m,uExa); // (Only for the extension of the project)
  
  // 4. Solve problem
  double tol = 1e-9; // (Currently useless)
  int maxit = 1000;
  jacobi(p.A, p.b, uNum, m, tol, maxit);
  //ConjGrad(p.A, p.b, uNum, m, tol, maxit);
  
  // 5. Compute error and export fields
  Vector uErr = uNum - uExa;
  saveToMsh(uNum, m, "solNum", "benchmark/solNum.msh");
  saveToMsh(uExa, m, "solRef", "benchmark/solExa.msh");
  saveToMsh(uErr, m, "solErr", "benchmark/solErr.msh");

  // Calcul de la norme L2 de l'erreur
  //Vector MuErr = p.M * uErr;
  //exchangeAddInterfMPI(MuErr, m);
  //Vector TuErr = uErr.transpose();
  //double L2Err = TuErr.dot(MuErr);
  // printf("\r Norme L2 = %e \n", L2Err);

  Vector MuExa = p.M * uExa;
  exchangeAddInterfMPI(MuExa, m);
  //Vector TuExa = uExa.transpose();
  //double L2Exa = TuExa.dot(MuExa);
  //double L2    = sqrt(L2Err/L2Exa);
  //printf("\r Norme L2 = %e \n", L2);

  Vector MuErr = p.M * uErr;
  exchangeAddInterfMPI(MuErr, m);
  
  /*double L2Err = 0; 
  for(int i=0; i<m.nbOfNodes; ++i){
    if(m.elemNoRep(i) > 0) {
      L2Err += MuErr(i)*uErr(i);
    }
  }
  double L2tot = L2Err;*/

  double L2tot = ParScalarProd(MuErr, uErr, m);
  double L2ex = ParScalarProd(MuExa, uExa, m);
  
  //double L2glob,L2globExa;
  //MPI_Reduce(&L2tot, &L2glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Reduce(&L2Exa, &L2globExa, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  double L2fin = sqrt(L2tot/L2ex);
  if(myRank == 0){
    printf("\r norme approchee L2 = %e \n", L2tot);
    printf("\r norme exa L2 = %e \n", L2ex);
    //printf("\r Norme L2 exa = %e \n", L2globExa);
    //printf("\r Norme L2 glob (ici somme des petites normes L2) = %e \n", L2glob);
    printf("\r Norme L2 final = %e \n", L2fin);
    }

  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
