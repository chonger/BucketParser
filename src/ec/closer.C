#include "CellClos.h"
#include "ECArgs.h"
#include "Tree.h"
#include <assert.h>
#include "Trellis.h"

static pthread_mutex_t readlock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t writelock = PTHREAD_MUTEX_INITIALIZER;
  
#define numThreads 1
float tcd=0;
float td=0;
float gd=0;

static void* mainloop(void* arg)
//void mainloop()
{
  int numSents=0;
  for(;;numSents++){
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
    int sz=yield.size();
    if(sz>=MAXLEN)continue;
    TAnses tAs;
    for(int i=0;i<=sz;i++){
      TAns ta;
      tAs.push_back(ta);
      
    }
    int pos=0;
    CellClos clClos;
    clClos.findTAnses(tree,tAs,pos);
    Sent st(yield,tAs);
    Trellis trel(st);

    clClos.closer(st);
    pthread_mutex_lock(&writelock);
    tcd+=clClos.totC;
    td+=clClos.tot;
    gd+=clClos.gld;
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
  CellClos::readFeats();
  //mainloop();

  int i;
  pthread_t thread[numThreads];
  for(i = 0 ; i < numThreads  ; i++){
    pthread_create(&thread[i],0,mainloop, 0);
  }
  for(i=0; i<numThreads; i++){
    pthread_join(thread[i],0);
  }

  float p=(tcd/td);
  float r=(tcd/gd);
  cerr<<p<<" "<<r<<" " <<(2*p*r)/(p+r)<<endl;
  exit(0);
}
