#include "CellClos.h"
#include "ECArgs.h"
#include "Tree.h"
#include <assert.h>
#include "Trellis.h"


void
goldLabs(Tree* t,Indicies& labs)
{
  Tree* p=t->subtrees;
  if(!p->subtrees){
    assert(t->label<OFFSET);
    labs.push_back(t->label);
    return;
  }
  for(;p;p=p->sibling)goldLabs(p,labs);
}


int
main(int argc, char *argv[])
{
  ECArgs args(argc, argv);
  SymbolTable ntst;
  SymbolTable tst;
  tst.offset=OFFSET;
  Tree::ntst=&ntst;
  Tree::tst=&tst;
  
  Sents sents;
  Trellis::init(ntst,tst);
  int numSents=0;
  CellClos clClos;
  clClos.doTraining=1;
  for(;;numSents++){
    Tree* tree=Tree::make(cin);
    if(!cin)break;
    if(!tree->label)continue;
    //if(numSents>10)break;
    tree->label=rootSym;
    Indicies gold;
    Indicies yield;
    tree->yield(yield);
    int sz=yield.size();
    for(int i=0;i<sz;i++) Trellis::words[yield[i]].cnt++;
    Trellis::words[stopSym].cnt++;
    if(sz>=MAXLEN)continue;
    goldLabs(tree,gold);

    TAnses tAs;
    for(int i=0;i<=sz;i++){
      TAns ta;
      tAs.push_back(ta);
      
    }
    int pos=0;
    clClos.findTAnses(tree,tAs,pos);
    Sent st(yield,tAs);
    Trellis trel(st);
    sents.push_back(st);
  }
  clClos.train(sents);
  clClos.writeFeats();
  exit(0);
}
