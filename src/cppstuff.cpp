#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <unistd.h>
#include <vector>
#include <queue>
#include <map>
#include <cstdio>
#include <limits>
#include "cppstuff.h"

using namespace std;

vector<vector<int> > graph;// the currently studied graph
vector<double> degreeDist; // degree distribution of the graph currently studied
vector<bool> isRed;        // isRed[i] is true iff sample i is red
vector<bool> isBlue;       // isBlue[i] is blue iff sample i is blue
int redTotal,blueTotal;    // total number or red (resp. blue) samples
int nGenes;                // total number of alterations
int nIndiv;                // total number of individuals
vector<int> globalHits;    // DO NOT TOUCH : functions step_score and make_children use this as a shared memory .

multiset<Node> waiting;    // nodes waiting to be added to the solution tree, in order of incerasing step-score
vector<Node> solutionTree; // structure containing the nodes of the solution tree
int leafCount=0;           // the number of leaves in the solution tree (is kept up to date during growth)
double currentThreshold=0; // the threshold is increased gradually during the growth of the solution tree
double goalThreshold=1;    // the threshold will not be increased past this value (can be set by user)
int sizeBound=10000;       // the solution tree will not grow past this size (can be set by user)

// For debugging : print a vector in this format : (23, 12, -3)
void printVector(const vector<int>& t)
{
  Rprintf("(");
  int i;
  for(i=0;i+1<t.size();i++) // CAUTION : size() is unsigned and could be zero, so do not compute t.size()-1 !
   Rprintf("%d, ",t[i]);
  if (t.size()>0)
   Rprintf("%d",t[i]);
  Rprintf(")");
}

//count occurrences of "true" in a vector
int boolWeight(const vector<bool> & t)
{
 int s=0;
 int l=t.size();
 for(int i=0;i<l;i++) if(t[i]) s++;
 return s;
}

// 'red' and 'blue' are the number of red/blue samples touched by some alteration set.
// 'hits' shows how many times each sample is hit by said alteration set
// 'gene' is an alteration we want to add :
// process this by modifying red, blue, and hits accordingly.
void processHits(int& red, int& blue, vector<int>& hits, int gene)
{
  for(vector<int>::iterator n=graph[gene].begin();n!=graph[gene].end();n++)
  {
   if(hits[*n]==0)
   { if(isRed[*n]) red++;
     if(isBlue[*n]) blue++;
   }
   hits[*n]++;
  }
}

// 'red' and 'blue' are the number of red/blue samples touched by some alteration set
// 'hits' shows how many times each sample is hit by said alteration set
// 'gene' is an alteration in this set, we want to REMOVE.
// process this by modifying red, blue, and hits accordingly.
void unProcessHits(int& red, int& blue, vector<int>& hits, int gene)
{
  for(vector<int>::iterator n=graph[gene].begin();n!=graph[gene].end();n++)
  {
   hits[*n]--;
   if(hits[*n]==0)
   { if(isRed[*n]) red--;
     if(isBlue[*n]) blue--;
   }
  }
}

// builds a score table (a scoreStruct object) designed to hold the values of score(r,b) for every r and every b
// where score(r,b) is the score of a solution hitting r red and b blue samples.
// AND also designed to hold score(g) for every gene g (computed only on demand, see individualScore(g) below).
scoreStruct::scoreStruct()
{
 memo=vector<vector<double> >();
 for(int r=0;r<=redTotal;r++)
 {memo.push_back(vector<double>());
  for(int b=0;b<=blueTotal;b++)
   memo[r].push_back(-R::phyper(b,blueTotal,redTotal,r+b, true,true)*r/redTotal); // (minus log) probability to draw at most b blue balls by drawing (r+b) balls at 
                                                                                  // random, TIMES the proportion of total red samples that we actually hit.
 }
 individualMemo=vector<double> (nGenes, -1);
}

// query the score table
double scoreStruct::operator()(int r, int b)
{ return memo[r][b];}

// get the score of a particular gene (by first counting its red/blue neighbors, then querying the score table)
// store the result in a memo of individual scores of single genes.
double scoreStruct::individualScore(int g)
{
 if(individualMemo[g]==-1)
 {
  int red=0,blue=0;
  for(vector<int>::iterator n=graph[g].begin();n!=graph[g].end();n++)
  { if(isRed[*n]) red++;
    if(isBlue[*n]) blue++;
  }
  individualMemo[g]=score(red,blue);
 }
 return individualMemo[g];
}

// distrStrict is a structure designed to contain many score *distributions*
// A score *distribution* is the repartition function of the score obtained by adding some random gene to a given alteration set
// The "current" distribution depends on how many red and blue samples are hit by the "current" set.
distrStruct::distrStruct()
{currentR=currentB=-1;}

// Query a distribution : get P(new score >= value) given that the current ses hits r red and b blue samples.
// If the current distribution has parameters r and b, just look it up.
// Otherwise, load the correct distribution (see loadDistr(r,b) ) to make it the current one, then look it up.
double distrStruct::operator()( int r, int b, double value)
{
 if (r != currentR || b != currentB) loadDistr(r,b);
 return (*scoresLaw).lower_bound(value)->second;
}

// load ONE of the score distributions you have in store :
// the one for adding to an alteration set hitting 'r' red and 'b' blue samples.
// If not in store yet, then "make" it : see makeDistr(r,b)
void distrStruct::loadDistr(int r, int b)
{
  pair<int,int> p=make_pair(r,b);
  map<pair<int,int>,multimap<double,double> >::iterator x=memo.find(p);
  if(x==memo.end()) x=makeDistr(p);
  scoresLaw=&(x->second);
  currentR=r;
  currentB=b;
}

// Compute one distribution : create it in the form of a multimap<double,double> and insert it with the key (r,b) in the big memo.
// It is a multimap because if there are more than one way to get some score value, this score value will appear more than once as a key.
// At first each score is inserted as a key, and the associated value is its probability.
// Then the multimap is iterated over backwards to accumulate the probabilities.
// in the end, the value associated with the first occurrence of key x, is the probability to get a score higher or equal to x
map<pair<int,int>,multimap<double,double> >::iterator distrStruct::makeDistr(pair<int,int>& p) 
{
  // TODO  make vectorized calls to R::dhyper, for faster computation ?
  // TODO save space by keeping in distr(r,b) only the scores higher than score(r,b) ? (not entirely sure this is safe)
  // TODO (more involved...) if memory space is short, keep track of how often each distribution is queried (this is much cheaper
  //      than the distribution itself) and only keep the most queried, while staying in the available memory.
  // create a map inside the map and get reference to it, as well as an iterator to be returned.
  map<pair<int,int>,multimap<double,double> >::iterator x=memo.insert(make_pair(p,multimap<double,double>())).first; // very simple :-)
  multimap<double,double>& law=x->second;
  // Save a few constants
  int r=p.first, b=p.second;
  int freeReds=redTotal-r;
  int freeBlues=blueTotal-b;
  int free=freeReds+freeBlues;
  // Now the law of the new samples count
  vector<double> newDegrees(free+1);
  for(int i=0;i<degreeDist.size();i++)
   if(degreeDist[i]>0.)
    for(int j=0;j<=i && j<=free;j++)
     newDegrees[j]+=degreeDist[i]*R::dhyper(j,free,nIndiv-free,i,false); // taking i samples among nindiv, get j free ones.
  // now the law of the new score
  for(int i=0;i<newDegrees.size() ;i++)
   if(newDegrees[i]>0.)
    for(int j=0;j<=freeReds && j<=i;j++)
     law.insert(make_pair(score(r+j,b+i-j),newDegrees[i]*R::dhyper(j,freeReds,freeBlues,i,false)));
  law.insert(make_pair(std::numeric_limits<double>::max(),0)); // Make sure law.lowerBound() always finds something.
  // now the cumulative law : law[value]=P(score>=value)
  double sum=0;
  for(multimap<double,double>::reverse_iterator i=law.rbegin();i!=law.rend();i++)
     sum=i->second=i->second+sum;
  /* // Print the distribution for debug purposes
  Rprintf("distribution partant de %d rouges, %d bleus \n",r,b);
  for(map<double,double>::reverse_iterator i=law.rbegin();i!=law.rend();i++)
     Rprintf("P( score >= %f )= %f \n",i->first,i->second);
  */
  return x;
}

//' Original Step-score
//'
//' Computes the original step-score from a solution to another.
//'
//' This function can only be called after \link{solutions} has been called
//' or at least \link{takeGraph}, to create the data structure which 
//' contain the score distributions. \cr
//' Only the number of red and blue samples altered in the start and finish
//' solutions are needed. The corresponding solutions do not have to 
//' actually exist.
//'
//' @param r0 The number of red samples to start from
//' @param b0 The number of blue samples to start from
//' @param r1 The number of red samples to get to
//' @param b1 The number of blue samples to get to
//' @return the step-score value. If the parameters have meaningless values 
//' (negative or greater than the number of red/blue samples) the result is -1.
//' @export
// [[Rcpp::export]]
float singleStepScore(int r0, int b0, int r1, int b1)
{
   if (r0<0 || b0<0 || r1<0 || b1<0 || r0>redTotal || b0>blueTotal || r1>redTotal || b1>blueTotal)
	   return -1;
   return distr(r0,b0,score(r1,b1));
}

//' Get Ready to Work with This Graph
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' Set Some Internal Global Variables to Deal with a Given Graph : 
//' the degree distribution of the alterations, the lists of neighbours,
//' the number of alterations, and the memo that will contain score
//' distribution tables used for step-score computation.
//'
//' @param data A vector of integer vectors. \code{data[[i]]} should be an integer vector containing the ids of the samples hit by alteration \code{i}
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void takeGraph(Rcpp::GenericVector & data)
{
 graph.clear();
 degreeDist.clear();
 nGenes=data.size();
 for(int i=0;i<nGenes;i++)
 {
  Rcpp::IntegerVector line=data[i];
  graph.push_back(Rcpp::as<vector<int> >(line));
  int l =line.size();
  while(degreeDist.size() <= l)
   degreeDist.push_back(0);
  degreeDist[l]++;
 }
 for(int i=0;i<degreeDist.size();i++)
   degreeDist[i]/=1.0*nGenes;
 distr=distrStruct(); // erase all previously saved distributions.
}

//' Get Ready to Work with These Groups
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' Set Some Internal Global Variables to Deal with a Given Grouping of
//' the Samples : number of red/blue samples and the groups themselves.
//'
//' @param red A logical vector. \code{red[[i]]} should be true iff sample \code{i} is red
//' @param blue A logical vector. \code{blue[[i]]} should be true iff sample \code{i} is blue
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void takeGroups(Rcpp::LogicalVector & red, Rcpp::LogicalVector & blue)
{
 nIndiv=red.size();
 globalHits=vector<int>(nIndiv);
 isRed=Rcpp::as<vector<bool> >(red);
 isBlue=Rcpp::as<vector<bool> >(blue);
 redTotal=boolWeight(isRed);
 blueTotal=boolWeight(isBlue);
 distr=distrStruct(); // erase all previously saved distributions.
 score=scoreStruct(); // erase all previously saved scores
}

// This one should be obvious
bool contains(const vector<int> & t, int x)
{
 for(vector<int>::const_iterator i=t.begin(); i!=t.end(); i++ ) 
   if (*i==x) return true;
 return false;
}

bool equals(const vector<int> & t, const vector<int>& u)
{
 if (t.size()!=u.size()) return false;
 for(vector<int>::const_iterator x=u.begin();x!=u.end(); x++)
  if( ! contains(t,*x) ) return false;
 return true;
}

// is another ordering of this same gene set already in the waiting queue or solution tree ?
// (invariant=true if every ordering of this gene set should have the same step-score. This eases the search for it)
bool alreadyFound(const Node& newNode, bool invariant)
{
 const vector<int> & target=newNode.genes;
/* printf("checking if ");
 printVector(target);
 printf(" has already been found or queued.\n"); */
 multiset<Node>::iterator begin= invariant ? waiting.lower_bound(newNode) : waiting.begin();
 multiset<Node>::iterator end= invariant ? waiting.upper_bound(newNode) : waiting.end();
 for(multiset<Node>::iterator n=begin; n!=end;n++)
  if ( equals( target, n->genes) )
   return true; // inclusion test succeeded : it is the same gene set.
 for(vector<Node>::iterator n=solutionTree.begin();n!=solutionTree.end();n++)
  if ( equals( target, n->genes) )
   return true; // inclusion test succeeded : it is the same gene set.
 /* else
  {
   printVector(target);
   printf(" != ");
   printVector(n-> genes);
   printf("\n");
  } */
 return false;
}

double stepScoreOriginal(const Node & parent, int newGene, Node & newNode)
{// TODO remove this line when confirmed useless (i.e. this function will be called only when the extension is better than the parent) 
 //if(score(newNode.nRed,newNode.nBlue)<=score(parent.nRed,parent.nBlue)) return 2; // forbid the extension if it's worse than the parent
 return distr(parent.nRed,parent.nBlue,score(newNode.nRed,newNode.nBlue));
}

// a first idea, but not recommended because its behavior is hard to define precisely when equivalent solutions have different step-scores : the first (chronologically) is kept.
double stepScoreOriginalUnique(const Node & parent, int newGene, Node & newNode)
{
 if(alreadyFound(newNode,false)) return 2;
 return stepScoreOriginal(parent,newGene,newNode);
}

// do not allow if newGene alone is better than some of parent's genes alone. (genes that are ex-aequo are compared by their id)
double stepScoreBestFirst(const Node & parent, int newGene, Node & newNode)
{
 double newGeneScore=score.individualScore(newGene);
 for(vector<int>::const_iterator g=parent.genes.begin();g!=parent.genes.end();g++)
  { double gscore=score.individualScore(*g);
    if(gscore<newGeneScore) return 2;
    if(gscore==newGeneScore && *g < newGene) return 2;// (genes that are ex-aequo are compared by their id)
  }
 return stepScoreOriginal(parent, newGene, newNode); // replaces the two lines below
 //if(score(newNode.nRed,newNode.nBlue)<=score(parent.nRed,parent.nBlue)) return 2;
 //return distr(parent.nRed,parent.nBlue,score(newNode.nRed,newNode.nBlue));
}

// look at the least significant step that could lead to this solution
double stepScoreStrict(const Node & parent, int newGene, Node & newNode)
{
 //first, check that the new gene is really new (this is implicitly enforced in the other step-score definitions, but not this one)
 if(contains(parent.genes, newGene)) return 2;
 // then remove the new gene to compute its step-score, 
 int red=newNode.nRed, blue=newNode.nBlue;
 double newScore=score(red,blue);
 unProcessHits(red,blue,globalHits,newGene);
 double maxStep=distr(red,blue,newScore); // remember the biggest 'step' value is the least significant.
 processHits(red,blue,globalHits,newGene);

 //then do the same with all the old ones, to find the worst step-score.
 bool badGene=false; // a bad gene is one such that this alteration set would be better without it. If we find one, just trash this alteration set.
 for(vector<int>::const_iterator g=parent.genes.begin();g!=parent.genes.end() && !badGene;g++)
 {
   unProcessHits(red,blue,globalHits,*g);
   if(newScore<=score(red,blue)) badGene=true; // check for bad gene
   double step=distr(red,blue,newScore);
   if(step>maxStep) maxStep=step;  // check for new worst step-score
   processHits(red,blue,globalHits,*g);
 }
 if(badGene) return 2;
 //now check for duplicate : to help this we already assign a stepscore to newNode
 newNode.step=maxStep;  // this stepscore will be a hint in the search for duplicates
 if(alreadyFound(newNode,true)) return 2;
 return maxStep;
}


double (*stepScore) (const Node & parent, int newGene, Node & newNode)=stepScoreStrict;

//construct a node to represent the empty set
Node::Node():nBlue(0),nRed(0),step(0),genes(),parentIndex(-1){}

// construct a node from a parent and a new gene. the red and blue neighbour counts are already computed by the caller who decides wether or not to create this node.
Node::Node(Node & parent, int newGene, int red, int blue):genes(parent.genes),parentIndex(-1),nRed(red),nBlue(blue)
{ genes.push_back(newGene);}

bool Node::operator<(const Node& other) const
{ return step<other.step;}

void Node::makeChildren(int index)// maybe this index should be a property of the parent node instead of an argument of this function...
{
 maxThreshold=1; // stepscore threshold above which this node will have real children.
 fill(globalHits.begin(),globalHits.end(),0); // start by gathering information about the start configuration
 int red=0,blue=0;
 for(vector<int>::iterator n=genes.begin();n!=genes.end();n++)
    processHits(red,blue,globalHits,*n);
 double parentScore=score(red,blue);
 for(int g=0;g<nGenes;g++)		// then see how to extend this configuration
 {
   if (g%1000==0) Rcpp::checkUserInterrupt();// do this at regular intervals
   //TODO early failure modes ? At least check that the new gene is not already in genes ! or better, check that the score is high enough to stand a chance in more involved stepmodes.
   // Note : this optimization could also reduce memory use for distributions, especially critical in STRICT mode.
   processHits(red,blue,globalHits,g); // check if.. (see below)
   Node newNode(*this,g,red,blue);
   double newStep=stepScore(*this,g,newNode); // the step score... 
   double newScore=score(red,blue); 		// and the new score...
   // ... are sufficient, i.e. the score is better, the stepscore is below the goalthreshold, and better than the end of the queue (or the queue is not full)
   if( newScore>parentScore && newStep<= goalThreshold && (waiting.size() < sizeBound || newStep<= waiting.rbegin()->step) )
   {			// ... and if they are, finalize and insert new node.
     newNode.step=newStep; // (stepscore might or might not be already set (depending on which score function was used). Anyway, set it now.)
     newNode.parentIndex=index;
     waiting.insert(newNode);
     maxThreshold=min(maxThreshold,newStep);
   }
   // remove the part of the queue that will never be popped.
   multiset<Node>::reverse_iterator last=waiting.rbegin();
   while (! waiting.empty() && waiting.size()-waiting.count(*last) >= sizeBound)     
   {  
      waiting.erase(*last);// erase MANY elements that have the same key
      last=waiting.rbegin();
   }
   unProcessHits(red,blue,globalHits,g);
 }
}

//' Perform One Growth Step of the Solution Tree
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' Increase the step-core threshold to the next possible value (which is the
//' step-score of the first node in the waiting queue). This allows one or 
//' many nodes to be added to the solution tree.
//'
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void growTree()
{
 if(waiting.empty()) 
 {
  Rprintf("Cannot grow : threshold reached.\n");
  return;
 }
 int index=solutionTree.size();
 currentThreshold=waiting.begin()->step;
 multiset<Node>::iterator next;
 while(!waiting.empty() && (next=waiting.begin())->step <= currentThreshold)
 {
  solutionTree.push_back(*next);
  waiting.erase(next);
  Node& node=solutionTree.back();
  node.minThreshold=currentThreshold;
  node.isLeaf=true;
  if(node.parentIndex!=-1 && solutionTree[node.parentIndex].isLeaf) // parent ceases to be a leaf, and child starts being one.
     solutionTree[node.parentIndex].isLeaf=false;
  else
     leafCount+=1;                                                  // if parent already wasn't a leaf, then we have one more.
  node.makeChildren(index++);
 } 
 Rprintf("threshold <- %f, Nodes : %ld \n",currentThreshold,solutionTree.size());
}


//' Create a Dataframe From the Solution Tree
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' Gather information about the nodes contained in it and write them to
//' A dataframe. Note that it only contains numbers, which will later
//' have to be translated by R code into names, unknown to C++ code.
//'
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
Rcpp::List encodeSolutionTree()
 {
  vector<vector<int> > families;
  vector<int> reds;
  vector<int> blues; 
  vector<int> parents; 
  vector<double> steps; 
  vector<double> scores; 
  vector<double> thresholds; 
  vector<double> maxThresholds; 
  vector<bool> leaf; 
  
  for(vector<Node>::iterator n=solutionTree.begin();n!=solutionTree.end();n++)
   {
    families.push_back(n->genes);
    reds.push_back(n->nRed);
    blues.push_back(n->nBlue);
    steps.push_back(n->step);
    scores.push_back(score(n->nRed,n->nBlue));
    thresholds.push_back(n->minThreshold);
    maxThresholds.push_back(n->maxThreshold);
    leaf.push_back(n->isLeaf);
    parents.push_back(n->parentIndex+1);
   }
  
  return Rcpp::List::create(Rcpp::Named("alterations",families),Rcpp::Named("reds",reds),Rcpp::Named("blues",blues),Rcpp::Named("score",scores),Rcpp::Named("step_score",steps),Rcpp::Named("threshold",thresholds),Rcpp::Named("parent",parents),Rcpp::Named("leaf",leaf),Rcpp::Named("childrenThreshold",maxThresholds));
 }

//' Pass a Step-score Mode to the C++ Code
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' @param stepMode : an integer. 0 for "original", 1 for "best-first", 2 for "strict"
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void takeStepScore(int stepMode)
{
 switch(stepMode)
 {case 0 : stepScore=stepScoreOriginal; break;
  case 1 : stepScore=stepScoreBestFirst; break;
  case 2 : stepScore=stepScoreStrict; break;
  case 3 : stepScore=stepScoreOriginalUnique; break; //discouraged
 }
}

//' Reset the Solution Tree
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' Create a fresh solution tree, empty and ready to be grown, 
//' clearing the previously built solution tree, if necessary.
//'
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void resetSolutionTree()
{
 solutionTree.clear();
 waiting.clear();
 waiting.insert(Node());
 currentThreshold=0;
 leafCount=0;
}

//' Pass a Tree Size Bound to the C++ code
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' @param treeSize An integer
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void setBound(int treeSize)
{ sizeBound=treeSize;}

//' Pass a Step-score Threshold to the C++ Code
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' @param threshold : An integer
//' @return nothing.
//' @family inner functions
//' @export
// [[Rcpp::export]]
void setThreshold(double threshold)
{ goalThreshold=threshold;}

//' Grow Solution Tree Until a Criterion is Met
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' The solution tree is grown using \code{\link{growTree}} until
//' its size reaches the set bound, or the step-score threshold reaches
//' the predefined value. 
//'
//' @return nothing.
//' @seealso \code{\link{growTree}}, \code{\link{setThreshold}}, \code{\link{setBound}}
//' @family inner functions
//' @export
// [[Rcpp::export]]
void fullyGrowTree()
{
 while (!waiting.empty() && solutionTree.size() < sizeBound)
  growTree();
 if (solutionTree.size()>= sizeBound) Rprintf("Size bound reached\n");
 if (waiting.empty()) Rprintf("Threshold reached\n");
}


//' Grow Solution Tree Until a Criterion is Met, forget the Solutions, keep the Growth history.
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' The solution tree is grown using \code{\link{growTree}} until
//' its size reaches the set bound, or the step-score threshold reaches
//' the predefined value. 
//'
//' @return A list showing the values successively taken by the step-score threshold, along
//' with the number of nodes and the number of leaves at each stage.
//' @seealso \code{\link{growTree}}, \code{\link{setThreshold}}, \code{\link{setBound}}
//' @family inner functions
//' @export
// [[Rcpp::export]]
Rcpp::List exploreTree()
{
  vector<int> nodes; 
  vector<int> leaves; 
  vector<double> thresholds; 

  resetSolutionTree();
  while (!waiting.empty() && solutionTree.size() < sizeBound)
  {
    growTree();
    thresholds.push_back(currentThreshold);
    nodes.push_back(solutionTree.size());
    leaves.push_back(leafCount);
  }
  return Rcpp::List::create(Rcpp::Named("threshold",thresholds),Rcpp::Named("nodes",nodes),Rcpp::Named("leaves",leaves));
}


//' Compute the Domination Relationship in the Current Graph
//'
//' (See details below about what "domination" means) 
//' Compute for every alteration A :
//'   \itemize{
//'     \item the list of alterations dominated by A "loosely" (i.e. including those dominating A)
//'     \item the list of alterations dominated by A "strictly" (i.e. excluding those dominating A)
//'     \item the list of alterations dominating A "loosely" (i.e. including those dominated by A)
//'     \item the list of alterations dominating A "strictly" (i.e. excluding those dominated by A)
//'   }
//' 
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//' 
//' Given a network connecting samples and alterations, and given two groups of samples, red and blue, alteration A
//' will dominate alteration B if : 
//' \itemize{
//'  \item At least <redPercent>\% of the red samples affected by B are also affected by A, 
//'  \item At least <bluePercent>\% of the blue samples affected by A are also affected by B,
//'  \item The loci of A and B lie on the same chromosome, at most <distance> base pairs apart (this part is for alterations which are genomic, and can simply be ignored for other applications)
//' }
//' This is a somewhat more flexible version of the more natual criteria according to which A would dominate B if : 
//' \itemize{
//'  \item The red samples affected by B are also affected by A
//'  \item The blue samples affected by A are also affected by B
//' } 
//' To compute this simpler and more restrictive relation, use \code{redPercent=100} and \code{bluePercent=100}.
//' 
//' @param redPercent The required percentage of red samples
//' @param bluePercent The required percentage of blue samples
//' @param chromosome (for genomic alterations) An integer vector containing the chromosome of each
//' alteration in the current graph (see \code{\link{takeGraph}}). If omitted, this parameter has no effect.
//' @param locus (for genomic alterations) An integer vector with the same length as \code{chromosome}, containing
//' the 'locus' of each alteration (e.g. in base pairs).  If omitted, this parameter has no effect.
//' @param distance The maximum distance accepted between two alterations' loci for
//' domination to take place. If omitted, this parameter has no effect.
//' @return A list \code{l} of 6 named lists of character vectors. If \code{A} is the name of an alteration, then :
//' \itemize{
//' \item \code{l$loose_dominators[[A]]} contains the names of alterations "loosely" dominating \code{A}.
//' \item \code{l$loose_dominated[[A]]} contains the names of alterations "loosely" dominated by \code{A}.
//' \item \code{l$strict_dominators[[A]]} contains the names of alterations "strictly" dominating \code{A}.
//' \item \code{l$strict_dominated[[A]]} contains the names of alterations "strictly" dominated by \code{A}.
//' } 
//' @seealso \code{\link{takeGraph}}, \code{\link{domination.graph}}
//' @family inner functions
//' @export
// [[Rcpp::export]]
Rcpp::List rawDominationGraphFineTuned(Rcpp::IntegerVector chromosome, Rcpp::IntegerVector locus, int distance, double redPercent, double bluePercent) 
{
 vector<vector<int> > strict_dominators(nGenes, vector<int>() );
 vector<vector<int> > strict_dominated(nGenes, vector<int>() );
 vector<vector<int> > loose_dominators(nGenes, vector<int>() );
 vector<vector<int> > loose_dominated(nGenes, vector<int>() );
 vector<vector<int> > equivalent(nGenes, vector<int>() );
 vector<bool> hasA(nIndiv,false);
 vector<bool> hasB(nIndiv,false);

 for(int a=0; a < nGenes; a++)
 {
 	 fill(hasA.begin(),hasA.end(),false);
 	 for(vector<int>::iterator n=graph[a].begin(); n != graph[a].end(); n++)
 	 	  hasA[*n]=true;

 	 for(int b=a+1; b < nGenes; b++)
 	 {
/*
bool debug=false;
if(debug) sleep(1);
if(debug) printf("debugging the pair %d %d\n",a,b);
if(debug) printf("chromosomes are %d and %d \n",chromosome[a],chromosome[b]);
if(debug) sleep(1);
*/
		 if(chromosome[a]!=chromosome[b]) continue;
                 if(Rcpp::IntegerVector::is_na(chromosome[a])) continue;
/*
if(debug) printf("same chromosome\n");
if(debug) sleep(1);
*/
		 if(abs(locus[a]-locus[b])>distance) continue;
/*
if(debug) printf("not too far apart\n");
if(debug) sleep(1);
*/
 	 	 fill(hasB.begin(),hasB.end(),false);
 	 	 for(vector<int>::iterator n=graph[b].begin(); n != graph[b].end(); n++)
 	 	 	 hasB[*n]=true;
/*
if(debug) printf("marked B's individuals\n");
if(debug) sleep(1);
*/
 	 	 int table[2][2][2]; // table[isRed][hasA][hasB]= number of individuals of given color having (or not) A, and/or B, accordingly.
                 for(int i=0;i<8;i++)
                     table[i/4%2][(i/2)%2][i%2]=0;

 	 	 for(int i=0;i<nIndiv;i++)
                 {
		   if (isRed[i])   table[1][int(hasA[i])][int(hasB[i])]++; 
                   if (isBlue[i])  table[0][int(hasA[i])][int(hasB[i])]++;// Those 2 lines could be one only if red and blue formed a partition.
		 }
/*
if(debug) printf("computed contingencies : \n red %d common, %d specific to a, %d specific to b\n  blue  %d common, %d specific to a, %d specific to b\n ",table[1][1][1],table[1][1][0],table[1][0][1],table[0][1][1],table[0][1][0],table[0][0][1]);
if(debug) sleep(1);
*/
		 bool adominatesb=(score.individualScore(a) >= score.individualScore(b) // a has a better score than b
                     && table[0][1][1]*100>=bluePercent*(table[0][1][1]+table[0][1][0]) // common blues >= bluePercent * blues of dominator
                     && table[1][1][1]*100>=redPercent* (table[1][1][1]+table[1][0][1]));  // common red >= redPercent * reds of dominated
                 bool bdominatesa=(score.individualScore(a) <= score.individualScore(b) // b has a better score than a
                     && table[0][1][1]*100>=bluePercent*(table[0][1][1]+table[0][0][1]) // common blues >= bluePercent * blues of dominator
                     && table[1][1][1]*100>=redPercent* (table[1][1][1]+table[1][1][0]));  // common red >= redPercent * reds of dominated

                 if (adominatesb)
                 {  loose_dominated[a].push_back(b);
                    loose_dominators[b].push_back(a);
                    if (bdominatesa)
                    {  loose_dominated[b].push_back(a);
                       loose_dominators[a].push_back(b);
                       equivalent[a].push_back(b);
                       equivalent[b].push_back(a);
/*
if(debug) printf("A and B equivalent \n");
if(debug) sleep(1);
*/
                    }
                    else
                    {
                       strict_dominated[a].push_back(b);
                       strict_dominators[b].push_back(a);
/*
if(debug) printf("A dominates B\n");
if(debug) sleep(1);
*/
                    }
                 }
                 else if (bdominatesa)
                 {  loose_dominated[b].push_back(a);
                    loose_dominators[a].push_back(b);
                    strict_dominated[b].push_back(a);
                    strict_dominators[a].push_back(b);
/*
if(debug) printf("B dominates A\n");
if(debug) sleep(1);
*/
                 }
 	 }
 }
 return Rcpp::List::create(Rcpp::Named("loose_dominators",loose_dominators),Rcpp::Named("loose_dominated",loose_dominated),Rcpp::Named("strict_dominators",strict_dominators),Rcpp::Named("strict_dominated",strict_dominated));
}



//' Compute Connected Components in a Graph.
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' @param neighbours An adjacency list in the form of a vector of vectors.
//' @return An integer vector \code{C} such that \code{C[[i]]} is a node in
//' the same component as \code{i}, and \code{C[[i]]=C[[j]]} iff \code{i} and
//' \code{j} are in the same component.
//' @family inner functions
//' @export
// [[Rcpp::export]]
Rcpp::GenericVector connectedComponents(Rcpp::GenericVector neighbours) 
{
 int n=neighbours.size();
 vector<bool> seen(n,false);
 vector<vector<int> > components=vector<vector<int> >();
 vector<int> stack=vector<int>();
 for(int i=0;i<n;i++)
   if(not seen[i])
   {
    vector<int> comp=vector<int>(); //the component
    seen[i]=true;
    stack.push_back(i);
    while (! stack.empty())
    {
      int j=stack[stack.size()-1];
      stack.pop_back();
      comp.push_back(j);
      vector<int> neighbours_of_j=Rcpp::as<vector<int> >(neighbours[j]);
      for(vector<int>::iterator k=neighbours_of_j.begin();k!=neighbours_of_j.end();k++)
         if(not seen[*k])
         {
           seen[*k]=true;
           stack.push_back(*k);
         }
    }
    components.push_back(comp);
   }
 return Rcpp::wrap(components);
}

//' Compute Strongly Connected Components in a Directed Graph.
//'
//' This is an inner function, called by other R functions in this
//' package, but it is left available for the expert user or developer.
//'
//' @details Tarjan's algorithm is used.
//'
//' @param neighbours An adjacency list in the form of a vector of vectors.
//' @return An integer vector \code{C} such that \code{C[[i]]} is a node in
//' the same component as \code{i}, and \code{C[[i]]=C[[j]]} iff \code{i} and
//' \code{j} are in the same component.
//' @family inner functions
//' @export
// [[Rcpp::export]]
Rcpp::GenericVector  scc(Rcpp::GenericVector neighbours) 
{  // TARJAN's algorithm of strongly connected components.
 int n=neighbours.size();
 vector<int> stack=vector<int>();
 // along the DFS, keep track of :
 vector<bool> added(n,false); // has this node been added to the tree or is it just waiting in the stack
 vector<int> level(n,-1); // level this node belongs to (might be updated, the last time is is pushed is the one that counts)
 vector<int> father(n,-1); // father (same as above, the last time this value is set is the one that counts)
 vector<int> lowlink(n,-1); // highest level reachable in one step from a descendant of this node. If the node can't reach above itself, it's the root of its SCC
 vector<int> scc(n,-1); // which strongly connected component is it in ?
 vector<int> root=vector<int>(); // root[k] will be the root of SCC number k.
 int current_scc=0; // will be incremented every time we meet the root (last node) of a SCC.

bool debug=false;
int debug_target=-1;

 //printf("initialized stuff, ready to start, n=%d\n",n); 
 for(int i=0;i<n;i++)
 {
   if(! added[i])
   {
    stack.push_back(i);
    level[i]=0;
    lowlink[i]=0;
    while (! stack.empty())
    {
      int j=stack.back();
if (j==debug_target) debug=true;
if(debug) printf("found %d on the stack, its level is %d , ",j,level[j]);
if(debug) printf("its father is %d \n",father[j]);
      if( ! added[j]) // first time it is found. Need to build descendance
      {
if(debug) printf("stacking its neighbours : ");
        added[j]=true;
        int l=level[j];
        vector<int> neighbours_of_j=Rcpp::as<vector<int> >(neighbours[j]);
        for(vector<int>::iterator k=neighbours_of_j.begin();k!=neighbours_of_j.end();k++)
           if(! added[*k]) // future descendant
           {
if(debug) printf("%d, ",*k);
             stack.push_back(*k);
             father[*k]=j;
             level[*k]=l+1;
             lowlink[*k]=l+1;
           }
           else if (scc[*k]==-1 && level[*k]<lowlink[j]) // ancestor
             lowlink[j]=level[*k];
           else if (father[*k]==-1) // left edge : if that node had no father, it has found one.
             father[*k]=j; // this father will help compute a normal form.
if(debug) printf("and its lowlink is now %d \n",lowlink[j]);
      }
      else if (scc[j]==-1) // second time it is found : descendance has been dealt with
      {
if(debug) printf("its successors are gone, its lowlink is %d. \n Put it in component %d ",lowlink[j],current_scc);
        stack.pop_back();
        scc[j]=current_scc;
        if (lowlink[j]>=level[j]) // it is a root
        { 
          root.push_back(j);
          current_scc++;
if(debug) printf("and switch to next component since it is a root.");
        }
        else if (father[j]>=0 && lowlink[j]<lowlink[father[j]])
           lowlink[father[j]]=lowlink[j];
if(debug && father[j]>=0) printf(". Its father's lowlink is %d \n", lowlink[father[j]]);

if (j==debug_target) debug=false;
      }
      else
      {
        stack.pop_back(); // has been pushed another time, and already dealt with.
if(debug) printf("already handled. Ignoring it.\n");
      }
    }
 
   }
 }
 //printf("main loop done : scc number %d just started and empty, ncomp=%d\n",current_scc,(int)root.size());
 // now wrap up everything, telling each guy the root of its component and preparing the result to be returned.
 vector<int> normal_form=vector<int>(root);
 int ncomp=root.size();
 int k;
 for(k=ncomp-1;k>=0;k--) // follow topologigal order
   if (father[root[k]]!=-1)
     normal_form[k]=normal_form[scc[father[root[k]]]];
 //printf("normal form done \n");
 vector<int> long_normal=vector<int>(n,-1);
 vector<int> long_root=vector<int>(n,-1);
 for(k=0;k<n;k++)
 {
   long_root[k]=root[scc[k]]; 
   long_normal[k]=normal_form[scc[k]];
 }
 //printf("long vectors done \n");
 vector< vector<int> > result=vector< vector<int> >();
 result.push_back(long_root);
 result.push_back(long_normal);
 //printf("result done");
 return Rcpp::wrap(result);
}

