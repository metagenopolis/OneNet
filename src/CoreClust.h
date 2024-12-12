#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

// // C++ header files
using namespace std;
using namespace Rcpp;
//
//
// #include <iostream>
// #include <iomanip>
// #include <fstream>
// #include <complex>
// #include <algorithm>
// #include <string>
// #include <limits>
//
//
// // C header files
// #include <stdio.h>
// #include <stdlib.h>
// #include <limits.h>
// #include <float.h>
// #include <math.h>
// #include <string.h>
// #include <iostream>     // std::cout
// #include <algorithm>    // std::sort
// #include <vector>       // std::vector


struct Edge
{
    unsigned int Node1;
    unsigned int Node2;
    double weight;
};


class CompactGraph{   //Class for undirected graphs
  private:

    int NodesNb;
    int EdgesNb;
    int MaxEdgesNb;
    std::vector<Edge> Graph;
    std::vector<unsigned int> PT_NodesNgb;  //nodes neighbors partition table
    std::vector<unsigned int> NodesNgb;           //nodes neighbors
    std::vector<int> Labels;
    std::vector<unsigned int> NodesNext;  //temporary variables for labels propagation
    std::vector<unsigned int> NodesCurr;  //temporary variables for labels propagation

  protected:

  public:

    /// Constructor
    CompactGraph(int MaxEdgesNb);

    /// Destructor
    ~CompactGraph();

    /// public functions

    //test function
    //arma::mat create_S() {
     // arma::mat S(2, 3);
      //S.fill(0);
      //return S;
    //};
  //
    //Read the graph in an ascci file representing the connectivity matrix.
    //-> The graph is undirected (only the terms above the diagonal will be read).
    //-> The file is supposed to only contain the matrix and no further information
    //-> The matrix has a size NbNodes*NbNodes
    //-> Only the edges with a weight higher than Thresh are loaded
  virtual void ReadInMatrix(arma::mat & M, float Thresh=0.);

  //   //Read a graph in a csv ascii file, where each row must have the following structure: [Node1] [Node2] [weight]
  //   //-> The graph is undirected.
  //   //-> The file is supposed to only contain the matrix and no further information
  //   //-> Node values MUST be integers between 0 and NbNodes-1
  //   //-> Only the edges with a weight higher than Thresh are loaded
  //  // virtual void ReadInList(char * FileName,float Thresh=0.);
  //
    //Save the graph as a list having the following structure: [Node1] [Node2] [weight]
  virtual arma::mat SaveInList();

    //Initiate the graph as the minimum/maximum spanning tree of the CompactGraph InputGraph
  virtual void MaximumSpanningTree(CompactGraph * InputGraph);
  //   virtual void MinimumSpanningTree(CompactGraph * InputGraph);
  //
  //
  // //Compute the core clusters of InputGraph -- greey version that works on the initial graph directely
  // //-> The size of the CORE-Clusters is in [CoreMinSize,2*CoreMinSize[
  // //-> Parameter $xi$ is a threshold on the lowest edge weight that will be scanned.
  // virtual void CoreClustering_Greedy(unsigned int CoreMinSize,double xi=0.5);
  //
  //Compute the core clusters of InputGraph -- standard version that works on a maximum spanning tree
  //-> The minimal size of the CORE-Clusters is CoreMinSize
  virtual void CoreClustering(unsigned int CoreMinSize);  //in the paper CoreMinSize=tau

  //sort a graph
  virtual void SortByIncreasingWeight();
  virtual void SortByDecreasingWeight();
  //
  // //compute the Partition table of the current graph
  virtual void ComputePartitionTable();

  //allocate memory for the labels and set them to 0
  // -> Type==0: all nodes have labels = 0
  // -> Type==1: node i has label i+1
  virtual void InitLabels(int InitType=0);

  // //Show the labels
  // virtual void ShowLabels();

  //Save the labels
  virtual arma::mat SaveLabels();
  //
  // //Label propagation in connected nodes with the same label.
  // //WARNING:
  // //    this->InitLabels, this->PT_NodesNgb and this->NodesNgb must be initiated and coherent with this->Graph.
   virtual unsigned int PropagateLabel(unsigned int RootNode, int LabelNew);
  //
  // //getters
   virtual int Get_NodesNb();
   virtual int Get_EdgesNb();
};


RCPP_MODULE(CompactGraph) {
  class_<CompactGraph>("CompactGraph")
  .constructor<int>()
  ;
}
