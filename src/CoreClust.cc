///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///
///Command line to compile:  g++ -o CoreClust CoreClust.cc -O -lm
///
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "CoreClust.h"


///-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-
/// PART. 1                                 CLASSES DEFINITION AND UTILITIES
///-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-


//constructor
CompactGraph::CompactGraph(int MaxEdgesNb){
  this->NodesNb=0;
  this->EdgesNb=0;
  this->MaxEdgesNb=MaxEdgesNb;
  this->Graph.resize(MaxEdgesNb);
  this->PT_NodesNgb.resize(0);
  this->NodesNgb.resize(0);
  this->Labels.resize(0);
  this->NodesNext.resize(0);
  this->NodesCurr.resize(0);
};

//destructor
CompactGraph::~CompactGraph(void){
  this->NodesNb=0;
  this->EdgesNb=0;
  this->MaxEdgesNb=0;
  this->Graph.resize(0);
  this->PT_NodesNgb.resize(0);
  this->NodesNgb.resize(0);
  this->Labels.resize(0);
  this->NodesNext.resize(0);
  this->NodesCurr.resize(0);
};







//Read the graph in an connectivity matrix.
//-> The graph is undirected (only the terms above the diagonal will be read).
//-> The matrix has a size NbNodes*NbNodes
//-> Only the edges with a weight higher than Thresh are loaded
void CompactGraph::ReadInMatrix(arma::mat & M, float Thresh){

  float tmpFl;
  int N1,N2;
  int NbNodes;

  //1) init
  NbNodes = M.n_cols;
  N1=0;
  N2=0;
  this->NodesNb=NbNodes;

  //2) read the graph  (values stricly above the diagonal)
  while(N1<NbNodes && N2<NbNodes){
    tmpFl = M(N1,N2);

    if ((N1>N2)&&(tmpFl>Thresh)){
      if (this->EdgesNb<this->MaxEdgesNb){

        this->Graph[this->EdgesNb].Node1=N1;
        this->Graph[this->EdgesNb].Node2=N2;
        this->Graph[this->EdgesNb].weight=tmpFl;
        this->EdgesNb++;
      }
      else{
        Rcerr << "More edges than expected are found in the matrix." << endl;
      }
    }
    N1++;
    if (N1>=this->NodesNb){
      N2++;
      N1=0;
    }
  }
  //3) generate the 'partition table' representation of the graph
  this->ComputePartitionTable();

  //4) allocate memory for the labels and set all of them to 0
  this->InitLabels();
}



//Only get the number of variables in the connexion matrix
int GetVariablesNumberInMatrix(arma::mat  & M){

  int NbNodes;

  NbNodes = M.n_cols;
  return NbNodes;
}


//Initiate the graph as the maximum spanning tree of the CompactGraph InputGraph
void CompactGraph::MaximumSpanningTree(CompactGraph * InputGraph){
  int CurrEdge;
  unsigned int Tst_Node1;
  unsigned int Tst_Node2;
  double Tst_weight;
  unsigned int NbPropagatedLabels;

  //1) Sort the edges of InputGraph by decreasing weights
  InputGraph->SortByDecreasingWeight();

  //2) init the labels of InputGraph with  label i for node i
  InputGraph->InitLabels(1);

  //3) Compute the maximum spanning tree
  this->NodesNb=InputGraph->Get_NodesNb();
  this->EdgesNb=0;

  for (CurrEdge=0;CurrEdge<InputGraph->Get_EdgesNb();CurrEdge++){

    if (CurrEdge>=this->MaxEdgesNb){
      Rcerr << "More edges than expected are found in the matrix." << endl;
    }
    else{
      Tst_Node1=InputGraph->Graph[CurrEdge].Node1;
      Tst_Node2=InputGraph->Graph[CurrEdge].Node2;
      Tst_weight=InputGraph->Graph[CurrEdge].weight;
      //cout << InputGraph->Labels[Tst_Node2];
      //cout << Tst_Node1;
    if (InputGraph->Labels[Tst_Node1]!=InputGraph->Labels[Tst_Node2]){ //the two nodes are in different rgions => merge the regions and keep the edge
      NbPropagatedLabels=InputGraph->PropagateLabel(Tst_Node2,InputGraph->Labels[Tst_Node1]);

      this->Graph[this->EdgesNb].Node1=Tst_Node1;
      this->Graph[this->EdgesNb].Node2=Tst_Node2;
      this->Graph[this->EdgesNb].weight=Tst_weight;
      this->EdgesNb++;
      }
    }
  }
  //4) generate the 'partition table' representation of the graph
  this->ComputePartitionTable();

  //5) allocate memory for the labels and set all of them to 0
  this->InitLabels();

}



// //functions to sort CompactGraphs
bool SortEdgesByWeights (const Edge &lhs, const Edge &rhs) { return lhs.weight < rhs.weight; }
bool SortEdgesByInvWeights (const Edge &lhs, const Edge &rhs) { return lhs.weight > rhs.weight; }

//Compute the core clusters of the graph -- standard version that works on a maximum spanning tree
//-> The minimal size of the CORE-Clusters is CoreMinSize
void CompactGraph::CoreClustering(unsigned int CoreMinSize){
  int CurrEdge;
  unsigned int Tst_Node1;
  unsigned int Tst_Node2;
  int CurrentMaxLabel,LabelCOREcluster,CurrLabel;
  double Tst_weight;
  unsigned int NbPropagatedLabels,TotalNbPropagatedLabels;
  int i,j;

  //1) Sort the graph edges by decreasing weights
  this->SortByIncreasingWeight();

  //2) init the graph labels with label 0 for all nodes
  this->InitLabels(0);
  CurrentMaxLabel=0;
  LabelCOREcluster=-1;

  //3) find the CORE-clusters
  for (CurrEdge=0;CurrEdge<this->EdgesNb;CurrEdge++){
    if (CurrEdge>=this->MaxEdgesNb){
      Rcerr << "More edges than expected are found in the matrix." << endl;
    }
    else{
      Tst_Node1=this->Graph[CurrEdge].Node1;
      Tst_Node2=this->Graph[CurrEdge].Node2;
      //cout << Tst_Node1 << " " <<  Tst_Node2<< ":" << endl;
      if ((this->Labels[Tst_Node1]==this->Labels[Tst_Node2])&&(this->Labels[Tst_Node2]>=0)){ //the two nodes are currently in the same region and are not frozen => cut

        //step 1: set the label CurrentMaxLabel in Tst_Node2
        CurrentMaxLabel++;
        CurrLabel=this->Labels[Tst_Node2];
        this->Labels[Tst_Node2] = CurrentMaxLabel; //change the label of Tst_Node2
        TotalNbPropagatedLabels=1;

        //step 2: propagate CurrentMaxLabel from Tst_Node2 neighbors (except Tst_Node1)
        for (j=this->PT_NodesNgb[Tst_Node2];j<this->PT_NodesNgb[Tst_Node2+1];j++)
          if ((this->NodesNgb[j]!=Tst_Node1)&&(this->PT_NodesNgb[this->NodesNgb[j]]>=0)&&(this->Labels[this->NodesNgb[j]]==CurrLabel)){
            NbPropagatedLabels=this->PropagateLabel(this->NodesNgb[j],CurrentMaxLabel);   //this->NodesNgb[j] is a Tst_Node2's neighbor
            TotalNbPropagatedLabels+=NbPropagatedLabels;
          }
          //cout << "S2 - TotalNbPropagatedLabels: " << " " <<  TotalNbPropagatedLabels << endl;

          //step 3: check whether the size of the Tst_Node2 side fits for a CORE-cluster: if yes freeze it
          if ((TotalNbPropagatedLabels>=CoreMinSize)&&(TotalNbPropagatedLabels<2*CoreMinSize)){ // a CORE-Cluster was identified -> freeze it
            NbPropagatedLabels=this->PropagateLabel(Tst_Node2,LabelCOREcluster);
            //cout << "(Edge " << CurrEdge << "/" << this->EdgesNb << ") Freeze N2" << " " <<  LabelCOREcluster <<  "  -> " << NbPropagatedLabels  << "nodes" << endl;
            LabelCOREcluster--;
          }

          //step 4: increment CurrentMaxLabel and measure the size of the Tst_Node1 side
          CurrentMaxLabel++;
          TotalNbPropagatedLabels=this->PropagateLabel(Tst_Node1,CurrentMaxLabel);
          //cout << "S1 - TotalNbPropagatedLabels: " << " " <<  TotalNbPropagatedLabels << endl;

          //step 5: check whether the size of the Tst_Node1 side fits for a CORE-cluster: if yes freeze it
          if ((TotalNbPropagatedLabels>=CoreMinSize)&&(TotalNbPropagatedLabels<2*CoreMinSize)){ // a CORE-Cluster was identified -> freeze it
            NbPropagatedLabels=this->PropagateLabel(Tst_Node1,LabelCOREcluster);
            //cout << "(Edge " << CurrEdge << "/" << this->EdgesNb << ") Freeze N1" << " " <<  LabelCOREcluster <<  "  -> " << NbPropagatedLabels  << "nodes" << endl;
            LabelCOREcluster--;
          }
      }
    }
  }

  //4) Post-process the labels
  for (i=0;i<this->NodesNb;i++){
    CurrLabel=this->Labels[i];
    if (CurrLabel>0)
      this->Labels[i]=0;
    else
      this->Labels[i]=-CurrLabel;
  }

  //5) outputs
  Rcout << endl;
  Rcout << "Found CORE-clusters: " << endl;
  for (j=1;j<-LabelCOREcluster;j++){
    Rcout << "CORE-cluster "  << j <<  ": ";
    for (i=0;i<this->NodesNb;i++)
      if (this->Labels[i]==j)
        Rcout << i << " ";
      Rcout << endl;
  }
}

void CompactGraph::SortByIncreasingWeight(){
  std::sort (this->Graph.begin(), this->Graph.begin()+this->EdgesNb, SortEdgesByWeights);
}

 void CompactGraph::SortByDecreasingWeight(){
   std::sort (this->Graph.begin(), this->Graph.begin()+this->EdgesNb, SortEdgesByInvWeights);
 }

//compute the partition table of the current graph
void CompactGraph::ComputePartitionTable(){
  int i;
  unsigned int tmpUInt1,tmpUInt2;
  unsigned int N1,N2;

  //1) define PT_NodesNgb

  this->PT_NodesNgb.resize(this->NodesNb+1);  //the value at this->PT_NodesNgb.resize(this->NodesNb) is used to know the nb of neighbors of node this->NodesNb-1

  for (i=0;i<this->NodesNb+1;i++)
    this->PT_NodesNgb[i]=0;

  for (i=0;i<this->EdgesNb;i++){
    N1=this->Graph[i].Node1;
    N2=this->Graph[i].Node2;
    this->PT_NodesNgb[N1]++;
    this->PT_NodesNgb[N2]++;
  }

  tmpUInt1=this->PT_NodesNgb[0];
  this->PT_NodesNgb[0]=0;
  for (i=1;i<this->NodesNb+1;i++){
    tmpUInt2=this->PT_NodesNgb[i];
    this->PT_NodesNgb[i]=tmpUInt1+this->PT_NodesNgb[i-1];
    tmpUInt1=tmpUInt2;
  }

  //for (i=0;i<this->NodesNb;i++)
  //   cout << i << " " << this->PT_NodesNgb[i] << endl;
  //cout << endl;

  //4.2) define memory of NodesNgb

  std::vector<unsigned int> tmpVec;
  tmpVec.resize(this->NodesNb);
  std::copy (this->PT_NodesNgb.begin(), this->PT_NodesNgb.begin()+this->NodesNb, tmpVec.begin() );

  this->NodesNgb.resize(this->PT_NodesNgb[this->NodesNb]);

  for (i=0;i<this->EdgesNb;i++){
    N1=this->Graph[i].Node1;
    N2=this->Graph[i].Node2;

    this->NodesNgb[tmpVec[N1]]=N2;
    tmpVec[N1]++;
    this->NodesNgb[tmpVec[N2]]=N1;
    tmpVec[N2]++;
  }

  //for (i=0;i<this->PT_NodesNgb[this->NodesNb-1]+tmpUInt1;i++)
  //   cout << this->NodesNgb[i] << endl;

}


//allocate memory for the labels and set them to 0
// -> Type==0: all nodes have labels = 0
// -> Type==1: nodei has label i+1
void CompactGraph::InitLabels(int InitType){
  int i;

  //allocate memory for the labels
  this->Labels.resize(this->NodesNb);

  //define initial labels
  if (InitType==0)
    for (i=0;i<this->NodesNb;i++)
      this->Labels[i]=0;

  if (InitType==1)
    for (i=0;i<this->NodesNb;i++)
      this->Labels[i]=i+1;

  //also allocate memory for the temporary vectors that will help propagating the labels
  this->NodesNext.resize(this->NodesNb);
  this->NodesCurr.resize(this->NodesNb);

}


//Label propagation in connected nodes with the same label.
//WARNING:
//    this->InitLabels, this->PT_NodesNgb and this->NodesNgb must be initiated and coherent with this->Graph.
unsigned int CompactGraph::PropagateLabel(unsigned int RootNode, int LabelNew){
  unsigned int NodesNext_Nb;
  unsigned int NodesCurr_Nb;
  int i,j,k;
  unsigned int currNode,ngbNode;
  int LabelOrig;
  int Test;
  unsigned int NbPropagatedLabels;


  //1) init
  if (this->Labels[RootNode]!=LabelNew){
    this->NodesNext[0]=RootNode;
    NodesNext_Nb=1;
    LabelOrig=this->Labels[RootNode];
  }
  else{
    NodesNext_Nb=0;
  }

  NbPropagatedLabels=0;

  //2) propagation
  while (NodesNext_Nb>0){
    //2.1) copy next in curr
    //cout << NodesNext_Nb  << ": ";
    for (i=0;i<NodesNext_Nb;i++){
      this->NodesCurr[i]=this->NodesNext[i];
      //cout << this->NodesCurr[i]  << " ";
    }
    //cout << endl;
    NodesCurr_Nb=NodesNext_Nb;
    NodesNext_Nb=0;
    //2.2) set the labels and search for neighbors
    for (i=0;i<NodesCurr_Nb;i++){
      currNode=this->NodesCurr[i];

      this->Labels[currNode]=LabelNew;
      NbPropagatedLabels++;
      for (j=this->PT_NodesNgb[currNode];j<this->PT_NodesNgb[currNode+1];j++){
        ngbNode=this->NodesNgb[j];
        if (this->Labels[ngbNode]==LabelOrig){
          Test=1; k=0;
          while ((Test==1)&&(k<NodesNext_Nb)){if (this->NodesNext[k]==ngbNode) Test=0; k++;} //to avoid having doubles
          for (k=0;k<NodesCurr_Nb;k++){if (this->NodesCurr[k]==ngbNode) Test=0;} //to avoid going through the same node several times
          if (Test==1){
            this->NodesNext[NodesNext_Nb]=ngbNode;
            NodesNext_Nb++;
          }
        }
      }

    }

  }
  return NbPropagatedLabels;
}
//
// //Show the labels
// void CompactGraph::ShowLabels(){
//   int i;
//
//   cout << "Labels: "  << endl;
//   for (i=0;i<this->NodesNb;i++)
//     cout << "Node "  << i << ": " << this->Labels[i] << endl;
// }
//
// //getters
int CompactGraph::Get_NodesNb(){return this->NodesNb;}
int CompactGraph::Get_EdgesNb(){return this->EdgesNb;}
//
//
// //UTILITIES
//
int string_to_int(string s) {
return atoi(s.c_str());
 }

//Save the graph as a list having the following structure: [Node1] [Node2] [weight]
//Only the edges with a weight higher than Thresh are written
arma::mat CompactGraph::SaveInList(){
  int i;

  arma::mat S(this->EdgesNb, 3);
  S.fill(0);

  //2) write data
  for (i=0;i<this->EdgesNb;i++){
    S(i,0)=this->Graph[i].Node1;
    S(i,1)=this->Graph[i].Node2;
    S(i,2)=this->Graph[i].weight;
  }
  return(S);
}

//Save the labels
arma::mat CompactGraph::SaveLabels(){
  int i;
  arma::mat V(this->NodesNb, 2);
  V.fill(0);

  //int m=this -> NodesNb;

  //2) write data
  for (i=0;i<this->NodesNb;i++){
    V(i,0)=i;
    V(i,1)=this->Labels[i];
  }
  //3) close file
  return(V);
}
//
// ///-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-
// /// PART. 2                                   MANAGEMENT FUNCTIONS
// ///-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-+*-
//
// [[Rcpp::export]]
Rcpp::List principalFunction(arma::mat &M, int & MinCoreSize){
  int NbVar;
  // int MinCoreSize = string_to_int(in_MinCoreSize);
  arma::mat V;
  arma::mat S;
  Rcout << "Minimal size of the CORE-clusters: " << MinCoreSize << endl;


  NbVar=GetVariablesNumberInMatrix(M);
  CompactGraph Graph(NbVar*NbVar);
  CompactGraph MSP(NbVar*NbVar);

  // //1) STANDARD VERSION

  Graph.ReadInMatrix(M);
  MSP.MaximumSpanningTree(&Graph);
  MSP.CoreClustering(MinCoreSize);
  V=MSP.SaveLabels();
  S=MSP.SaveInList();

  return Rcpp::List::create(S,V);
}




