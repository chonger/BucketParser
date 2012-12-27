#include "ECArgs.h"
#include "ECString.h"
#include "Tree.h"
#include <assert.h>
#include <ctype.h>
#include "Trellis.h"

static pthread_mutex_t readlock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t writelock = PTHREAD_MUTEX_INITIALIZER;

#define numThreads 1

int totTags=0;
int totCorrect=0;

void
goldLabs(Tree* t,Indicies& labs)
{
  Tree* p=t->subtrees;
  if(!p->subtrees){
    labs.push_back(t->label);
    return;
  }
  for(;p;p=p->sibling)goldLabs(p,labs);
}

int
score(Indicies& gold, Indicies& ans)
{
  int cor=0;
  IndiciesIter iIa=ans.begin();
  for(IndiciesIter iIg=gold.begin();iIg!=gold.end();iIg++){
    if(*iIa==*iIg)cor++;
    iIa++;
  }
  assert(iIa==ans.end());
  return cor;
}

int
score2(Indicies& gold, IndiciesV& ans)
{
  //cerr<<"S2-------------"<<endl;
  int cor=0;
  IndiciesVIter iIa=ans.begin();
  int loc=0; 
  for(IndiciesIter iIg=gold.begin();iIg!=gold.end();iIg++){
    Indicies& ids=*iIa;
    int sz=ids.size();
    //cerr<<loc<<" "<<*iIg<<" "<<sz<<" ";
    for(int i=0;i<sz;i++){
      if(ids[i]==*iIg){
	cor++;
	break;
      }
      //cerr<<ids[i]<<" ";
    }
    //cerr<<endl;
    iIa++;
    loc++;
  }
  assert(iIa==ans.end());
  assert(cor<=gold.size());
  return cor;
}


static void* mainloop(void* arg)
{
  for(int ct=0;;ct++){
    pthread_mutex_lock(&readlock);
    Tree* tree=Tree::make(cin);
    if(!cin){
      pthread_mutex_unlock(&readlock);
      break;
    }
    pthread_mutex_unlock(&readlock);
    if(!tree->label)continue;
    tree->label=rootSym;
    Indicies gold;
    Indicies yield;
    tree->yield(yield);
    if(yield.size()>=MAXLEN)continue;
    goldLabs(tree, gold);
    Sent snt(yield);

    Trellis trel(snt);

    int tc=score2(gold,trel.multiresult);
    
    pthread_mutex_lock(&writelock);
    totTags+=trel.sz;
    totCorrect+=tc;
    pthread_mutex_unlock(&writelock);

  }
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
  Trellis::init(ntst,tst);
  pthread_t thread[numThreads];

  int i;
  for(i = 0 ; i < numThreads  ; i++){
    pthread_create(&thread[i],0,mainloop, 0);
  }
  for(i=0; i<numThreads; i++){
    pthread_join(thread[i],0);
  }
  cerr<<totCorrect<<" "<<totTags<<" "
    <<(double)totCorrect/(double)totTags<<endl;
  exit(0);
}
