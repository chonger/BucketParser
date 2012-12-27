
#ifndef TREE_H
#define TREE_H

#include "ECString.h"
#include <fstream>
#include <sys/resource.h>
#include <iostream>
#include <stdlib.h>
#include "SymbolTable.h"

#define NEWSTORY 9999999
#define CATSEP		"-=~"		/* separates categories */
#define MAXLABELLEN 512
typedef vector<int> Ruleis;

class Tree
{
 public:
  Tree(stIndex i) : label(i),head(NULL),pprob(0),alpha(0),subtrees(NULL),sibling(NULL),parent(NULL),pred(NULL)
    { }
  Tree(Tree& other);
  ~Tree(){
    if(subtrees) delete subtrees;
    if(sibling) delete sibling;
  }
  static void init(SymbolTable* ntst1,SymbolTable* tst1)
    {
      ntst=ntst1;
      tst=tst1;
    }
  static void readStory(istream& is);
  void yield(vector<stIndex>& sti);
  void writeYield(ostream& os);
  stIndex     label;
  Tree*     head;
  double       pprob;
  double       alpha;
  stIndex      ptHack;
  Tree  *subtrees;
  Tree  *sibling;
  Tree  *parent;
  Tree  *pred;
  Indicies fnVals;
  void  write(ostream& os);
  void  w(ostream& os);
  static Tree* make(istream& is);
  static int singleLineReader; //default=1, set to 0 to read prettyprinted trees.
  Tree* hTreeFromTree();
  Tree* rhs1(){return subtrees;}
  Tree* rhs2(){return subtrees->sibling;}
  static Tree* readB(istream& is,int& cp,string& buf);
  double beta(stIndex subi);
  static SymbolTable* ntst;
  static SymbolTable* tst;
  int  wPos();
  int wPos(Tree* ht,int& pos);
  static int toLower;
 private:
  void   writeVals(ostream& os);
};
 
void addParent(Tree* tr,Tree* par);
typedef vector<Tree*> Trees;
typedef Trees::iterator TreesIter;
typedef vector<int> Ints;
typedef Ints::iterator IntsIter;

#endif
