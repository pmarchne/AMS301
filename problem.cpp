#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Compute the matrices of the linear system
//================================================================================

void buildLinearSystem(Problem& p, Mesh& m, double alpha, Vector& f)
{
  if(myRank == 0)
    printf("== build linear system\n");
  
  p.K.resize(m.nbOfNodes, m.nbOfNodes);
  p.M.resize(m.nbOfNodes, m.nbOfNodes);
  
  for(int iTri=0; iTri<m.nbOfTri; ++iTri){
    
    IntVector s = m.triNodes.row(iTri);
    Vector s0 = m.coords.row(s(0));
    Vector s1 = m.coords.row(s(1));
    Vector s2 = m.coords.row(s(2));
    
    double x01 = s0(0)-s1(0);
    double y01 = s0(1)-s1(1);
    double x12 = s1(0)-s2(0);
    double y12 = s1(1)-s2(1);
    double x20 = s2(0)-s0(0);
    double y20 = s2(1)-s0(1);
    
    double a = 0.5*abs(x12*y20 - x20*y12);
    
    double mel[3][3];
    mel[0][0] = a/6.;
    mel[0][1] = a/12.;
    mel[0][2] = a/12.;
    mel[1][0] = a/12.;
    mel[1][1] = a/6.;
    mel[1][2] = a/12.;
    mel[2][0] = a/12.;
    mel[2][1] = a/12.;
    mel[2][2] = a/6.;
    
    double kel[3][3];
    kel[0][0] = (x12*x12+y12*y12)/(4.*a);
    kel[0][1] = (x12*x20+y12*y20)/(4.*a);
    kel[0][2] = (x12*x01+y12*y01)/(4.*a);
    kel[1][0] = (x12*x20+y12*y20)/(4.*a);
    kel[1][1] = (x20*x20+y20*y20)/(4.*a);
    kel[1][2] = (x20*x01+y20*y01)/(4.*a);
    kel[2][0] = (x12*x01+y12*y01)/(4.*a);
    kel[2][1] = (x20*x01+y20*y01)/(4.*a);
    kel[2][2] = (x01*x01+y01*y01)/(4.*a);
    
    for(int j=0; j<3; ++j){
      for(int k=0; k<3; ++k){
        p.K.coeffRef(s(j), s(k)) += kel[j][k];
        p.M.coeffRef(s(j), s(k)) += mel[j][k];
      }
    }
  }
  
  p.A = alpha*p.M + p.K;
  p.b = p.M*f;
  exchangeAddInterfMPI(p.b, m);
}

//================================================================================
// Special treatment of the matrices for a Dirichlet boundary condition
//================================================================================

void buildDirichletBC(Problem& p, Mesh& m, Vector& uExa)
{
  if(myRank == 0)
    printf("== build Dirichlet BC\n");
  
  // Build 'relevement' function and mask (0 = interior node, 1 = boundary node)
  Vector g(m.nbOfNodes);
  IntVector isBnd(m.nbOfNodes);
  for(int i=0; i<m.nbOfNodes; ++i){
    g(i) = 0.;
    isBnd(i) = 0;
  }
  for(int iLin=0; iLin<m.nbOfLin; ++iLin){
    int i0 = m.linNodes(iLin,0);
    int i1 = m.linNodes(iLin,1);
    g(i0) = uExa(i0);
    g(i1) = uExa(i1);
    isBnd(i0) = 1;
    isBnd(i1) = 1;
  }
  
  // Modification of the system with 'relevement'
  Vector Ag = p.A*g;
  exchangeAddInterfMPI(Ag, m);
  p.b = p.b - Ag;
  
  // Pseudo-reduction
  for(int iTri=0; iTri<m.nbOfTri; ++iTri){
    IntVector s = m.triNodes.row(iTri);
    for(int j=0; j<3; ++j){
      for(int k=0; k<3; ++k){
        if((isBnd(s(j)) == 1) || (isBnd(s(k)) == 1)){
          if(j == k){
            p.A.coeffRef(s(j), s(k)) = 1.;
            p.b(s(j)) = g(s(j));
          }
          else {
            p.A.coeffRef(s(j), s(k)) = 0.;
          }
        }
      }
    }
  }
}
