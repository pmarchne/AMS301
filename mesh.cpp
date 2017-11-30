#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Read the mesh from a gmsh-file (.msh) and store in a mesh-structure 'm'
//================================================================================

void readMsh(Mesh& m, string fileName)
{
  if(myRank == 0)
    printf("== read mesh file\n");
  
  // Open mesh file
  ifstream meshfile(fileName.c_str());
  if(!meshfile){
    printf("ERROR: Mesh '%s' not opened.\n", fileName.c_str());
    exit(EXIT_FAILURE);
  }
  
  // Read mesh file
  double dummy;
  while(meshfile.good()){
    
    string line;
    getline(meshfile,line);
    
    // Read NODES
    if(line.compare("$Nodes") == 0){
      
      // Read the number of nodes
      meshfile >> m.nbOfNodes;
      
      // Read the node coordinates
      m.coords.resize(m.nbOfNodes,3);
      m.elemNoRep.resize(m.nbOfNodes);
      for(int i=0; i<m.nbOfNodes; ++i)
        meshfile >> dummy >> m.coords(i,0) >> m.coords(i,1) >> m.coords(i,2);
    }
    
    // Read ELEMENTS
    if(line.compare("$Elements") == 0){
      
      // Read the total number of elements
      int nbOfElems;
      meshfile >> nbOfElems;
      
      // Temporary arrays for elements infos/nodes
      IntVector elemNum(nbOfElems);          // gmsh numerbing
      IntVector elemType(nbOfElems);         // type of element
      IntMatrix elemNodes(nbOfElems,3);      // associated nodes (max 3. for triangle)
      IntVector elemPart(nbOfElems);         // partition tag
      
      // Read infos/nodes for all the elements
      m.nbOfLin = 0;
      m.nbOfTri = 0;
      for(int i=0; i<nbOfElems; i++){
        
        getline(meshfile, line);
        int infos;
        
        // Save element infos
        meshfile >> elemNum(i);
        meshfile >> elemType(i);
        meshfile >> infos;
        meshfile >> dummy >> dummy;
        if(infos > 2)
          meshfile >> dummy >> elemPart(i); // partition tag
        else
          elemPart(i) = 0;
        for(int j=5; j<=infos; j++)
          meshfile >> dummy;                // useless infos
        
        // Check element type and save element nodes
        switch (elemType(i)){
          case 1:  // segment
            m.nbOfLin++;
            meshfile >> elemNodes(i,0);
            meshfile >> elemNodes(i,1);
            break;
          case 2:  // triangle
            m.nbOfTri++;
            meshfile >> elemNodes(i,0);
            meshfile >> elemNodes(i,1);
            meshfile >> elemNodes(i,2);
            break;
          default:
            printf("ERROR: Element type '%i' not supported.\n", elemType(i));
            exit(EXIT_FAILURE);
            break;
        }
      }
      
      // Resize arrays for elements infos/nodes for each type of element
      if(m.nbOfLin > 0){
        m.linNum.resize(m.nbOfLin);
        m.linNodes.resize(m.nbOfLin,2);
        m.linPart.resize(m.nbOfLin);
      }
      if(m.nbOfTri > 0){
        m.triNum.resize(m.nbOfTri);
        m.triNodes.resize(m.nbOfTri,3);
        m.triPart.resize(m.nbOfTri);
      }
      
      // Save element infos/nodes for each type of element
      int iLin = 0;
      int iTri = 0;
      for(int i=0; i<nbOfElems; i++){
        switch (elemType(i)){
          case 1:  // line
            m.linNum(iLin) = elemNum(i);
            m.linNodes(iLin,0) = elemNodes(i,0)-1; // (gmsh is 1-index, here is 0-index)
            m.linNodes(iLin,1) = elemNodes(i,1)-1;
            m.linPart(iLin) = (elemPart(i)-1) % nbTasks; // (gmsh is 1-index, here is 0-index)
            iLin++;
            break;
          case 2:  // triangle
            m.triNum(iTri) = elemNum(i);
            m.triNodes(iTri,0) = elemNodes(i,0)-1;
            m.triNodes(iTri,1) = elemNodes(i,1)-1;
            m.triNodes(iTri,2) = elemNodes(i,2)-1;
            m.triPart(iTri) = (elemPart(i)-1) % nbTasks; // (gmsh is 1-index, here is 0-index)
            iTri++;
            break;
          default:
            printf("ERROR: Element type '%i' not supported.\n", elemType(i));
            exit(EXIT_FAILURE);
            break;
        }
      }
    }
  }
  
  meshfile.close();
  
  if(myRank == 0){
    printf("   -> %i nodes\n", m.nbOfNodes);
    printf("   -> %i boundary lines\n", m.nbOfLin);
    printf("   -> %i triangles\n", m.nbOfTri);
  }
}

//================================================================================
// Save a solution 'u' in a mesh file (.msh)
//================================================================================

void saveToMsh(Vector& u, Mesh& m, string name, string fileName)
{
  if(nbTasks > 1){
    ostringstream ss;
    ss << fileName << "_" << myRank;
    fileName = ss.str();
  }

  ofstream posFile(fileName.c_str());
  posFile << "$MeshFormat" << endl;
  posFile << "2.1 0 8" << endl;
  posFile << "$EndMeshFormat" << endl;
  
  posFile << "$ElementNodeData" << endl;
  posFile << "2" << endl;
  posFile << "\"" << name.c_str() << "\"" << endl;  // name of the view
  posFile << "" << endl;
  posFile << "1" << endl;
  posFile << "0" << endl; // ("Time")
  if(nbTasks > 0)
    posFile << "4" << endl;
  else
    posFile << "4" << endl;
  posFile << "0" << endl; // ("timeStep")
  posFile << "1" << endl; // ("numComp")
  posFile << m.nbOfTri << endl;   // total number of elementNodeData in this file
  if(nbTasks > 0)
    posFile << myRank << endl;
  for(int iTri=0; iTri<m.nbOfTri; iTri++){
    posFile << m.triNum[iTri] << " " << 3;
    for(int v=0; v<3; v++){
      int i = m.triNodes(iTri,v);
      posFile << " " << u(i);
    }
    posFile << endl;
  }
  posFile << "$EndElementNodeData" << endl;
  posFile.close();
}
