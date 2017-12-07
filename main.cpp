#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  // 0. Input arguments
  double alpha   = atof(argv[1]);  
  int    BC      = atoi(argv[2]); // 0 for Neumann BC and 1 for Dirichlet BC (Boundary Condition)
  int    Solver  = atoi(argv[3]); // 0 for Jacobi and 1 for CG (Conjugate Gradient)
  int    maxit   = atoi(argv[4]);
  //int    mesh_nb = atoi(argv[5]);

  // 1. Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
  
  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  Mesh m;
  readMsh(m, "benchmark/carre_h0.025_N4.msh");
  buildListsNodesMPI(m);
  buildLocalNumbering(m);
  
  // 3. Build problem (fields and system)
  Vector uNum(m.nbOfNodes);
  Vector uExa(m.nbOfNodes);
  Vector f(m.nbOfNodes);
  //double alpha = 1;
  for(int i=0; i<m.nbOfNodes; ++i){
    double x = m.coords(i,0);
    double y = m.coords(i,1);
    uNum(i) = 0.;
    switch(BC)
    {
      case 0 :
        uExa(i) = cos(M_PI*x)*cos(2*M_PI*y); // Neumann
        break;
      case 1 :
	uExa(i) = sin(M_PI*x)*sin(2*M_PI*y); // Dirichlet
	break;
    }
    f(i) = (alpha+5*pow(M_PI,2)) * uExa(i);
  }
  
  Problem p;
  buildLinearSystem(p,m,alpha,f);
  if (BC == 1)
    buildDirichletBC(p,m,uExa); // If Dirichlet
  
  // 4. Solve problem
  double tol = 1e-9; // (Currently useless)
  switch(Solver)
  {
    case 0 :
      jacobi(p.A, p.b, uNum, m, tol, maxit);   // Jacobi
      break;
    case 1 :
      ConjGrad(p.A, p.b, uNum, m, tol, maxit); // Conjugate Gradient
      break;
  }
  
  // 5. Compute error and export fields
  Vector uErr = uNum - uExa;
  saveToMsh(uNum, m, "solNum", "benchmark/solNum.msh");
  saveToMsh(uExa, m, "solRef", "benchmark/solExa.msh");
  saveToMsh(uErr, m, "solErr", "benchmark/solErr.msh");

  Vector MuExa = p.M * uExa;
  exchangeAddInterfMPI(MuExa, m);

  Vector MuErr = p.M * uErr;
  exchangeAddInterfMPI(MuErr, m);

  double L2tot = ParScalarProd(MuErr, uErr, m);
  double L2ex = ParScalarProd(MuExa, uExa, m);
  
  double L2fin = sqrt(L2tot/L2ex);
  if(myRank == 0){
    /*printf("\r norme approchee L2 = %e \n", L2tot);
    printf("\r norme exa L2 = %e \n", L2ex);*/
    printf("\r Norme L2 final = %e \n", L2fin);
    }

  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
