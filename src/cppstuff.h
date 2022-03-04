#include <vector>

using namespace std;

struct scoreStruct
{
  vector<vector<double> > memo;
  vector<double>  individualMemo;
  scoreStruct();
  double operator()(int,int);
  double individualScore(int);
};

struct distrStruct
{ int currentR, currentB;
  map< pair<int,int>, multimap<double,double> > memo;
  multimap<double,double>* scoresLaw;
  distrStruct();
  double operator()(int,int, double);
  void loadDistr(int,int);
  map<pair<int,int>,multimap<double,double> >::iterator makeDistr(pair<int,int>&);
};

scoreStruct score;
distrStruct distr;


struct PairComp
{ // best solutions are kept at the end : they have higher scores
  bool operator()(pair<int,int> a,pair<int,int>b) 
  {return score(a.first,a.second)<score(b.first,b.second);
  }
} comp;

vector< multimap< pair<int,int>,vector<int>, PairComp > > solutions;

struct Node
{
 vector<int> genes;
 int nBlue;
 int nRed;
 double step;
 int parentIndex;
 bool isLeaf;
 double minThreshold; // threshold must be below this for the node to exist
 double maxThreshold; // if threshold is above this, the node has children
 Node();
 Node             (Node& parent, int newGene, int red, int blue);// constructor
// void reinitialize(Node& parent, int newGene, int red, int blue);// equivalent but from already allocated object.
 bool operator<(const Node& other)const;
 void makeChildren(int index);
};

