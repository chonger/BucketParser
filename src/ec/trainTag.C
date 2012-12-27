#include "ECArgs.h"
#include "ECString.h"
#include "Tree.h"
#include <assert.h>
#include <ctype.h>
#include "Unkify.h"

#define NUMSTATES 49
#define MAXNUMWORDS 100000
stIndex rootSym;
stIndex stopSym;
stIndex stopWord;
int ntrees=0;
int wordLimit=0;

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

class Word
{
public:
  Word():cnt(0){for(int i=0;i<50;i++)tags[i]=0;}
  int cnt;
  int tags[50];
};

class State
{
public:
  State()
  {for(int i=0;i<50;i++){ntags[i]=0;for(int j=0;j<50;j++)nntags[i][j]=0;}}
  int ntags[50];
  int nntags[50][50];
};

Word words[MAXNUMWORDS];
State states[50];

void
printStats()
{
  for(int i=1;i<NUMSTATES;i++){
    State& st=states[i];
    for(int j=1;j<NUMSTATES;j++){
      for(int k=1;k<NUMSTATES;k++){
	cout<< Tree::ntst->toString(i)<<" "<<
	  Tree::ntst->toString(j)<<" "<<
	  Tree::ntst->toString(k)<<" "<<
	  st.nntags[j][k]<<endl;
      }
    }
  }
  for(int i=OFFSET+1;i<wordLimit;i++){
    Word& w=words[i];
    int j=0;
    int sum=0;
    for(;j<NUMSTATES;j++){sum+=w.tags[j];}
    if(sum==0&&!isUnk(i))continue;
    cout<<Tree::tst->toString(i)<<"\t";
    for(int j=0;j<50;j++){if(w.tags[j]>0)cout<<j<<" "<<w.tags[j]<<" ";}
    cout<<9999999<<endl;
  }
}

void
countWords(Tree* t)
{
  if(!t->subtrees)    words[t->label].cnt++;
  else for(Tree* p=t->subtrees;p;p=p->sibling)countWords(p);
}


void
binarize(Tree* t,stIndex& tgm2,stIndex& tgm1)
{
  Tree* p=t->subtrees;
  if(p->subtrees){
    for(;p;p=p->sibling)    binarize(p, tgm2,tgm1);
    return;
  }
  if(words[p->label].cnt<=3){
    p->label=unkify(p->label);
  }
  Word& w=words[p->label];
  w.tags[t->label]++;
  State& st=states[tgm2];
  st.ntags[tgm1]++;
  st.nntags[tgm1][t->label]++;
  stIndex xx1=tgm1;
  tgm2=xx1;
  tgm1=t->label;
}


void
readSyms()
{
  ifstream ntfs(Parser::dataDir + "nTerms.txt");
  for(;;){
    string nt;
    ntfs>>nt;
    if(!ntfs)break;
    Tree::ntst->toIndex(nt);
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
  readSyms();
  rootSym=ntst.toIndex("Root");
  stopSym=ntst.toIndex("STOP");
  stopWord=tst.toIndex("sToP");
  
  Trees trees;
  for(;;){
    Tree* tree=Tree::make(cin);
    if(!cin)break;
    if(!tree->label)continue;
    tree->label=rootSym;
    trees.push_back(tree);
    countWords(tree);
    words[stopWord].tags[stopSym]++;
  }
  for(TreesIter tI=trees.begin();tI!=trees.end();tI++){
    Tree* t=*(tI);
    stIndex nt1=stopSym;
    stIndex nt2=stopSym;
    binarize(t,nt1,nt2);
    State& st=states[nt1];
    st.ntags[nt2]++;
    st.nntags[nt2][stopSym]++;
    State& st2=states[nt2];
    st2.ntags[stopSym]++;
    st2.nntags[stopSym][stopSym]++;
    ntrees++;
    //cerr<<"TSTSZ "<<ntst.size()<<" "<<ntst.toString(ntst.size()-2)<<endl;
  }
  wordLimit=tst.size()+OFFSET;
  printStats();
  exit(0);
}
/*
stIndex
unkify(stIndex w)
{
  string ws=Tree::tst->toString(w);
  string uk="UNK";
  int sz=ws.size()-1;
  int i=0;
  if(isupper(ws[0]))uk="C"+uk;
  if(isdigit(ws[0])&&isdigit(ws[sz]))uk=uk+"N";
  else if(sz<=2){}
  else if(ws[sz]=='g'&&ws[sz-1]=='n'&&ws[sz-2]=='i')uk=uk+"ING";
  else if(ws[sz]=='d'&&ws[sz-1]=='e')uk=uk+"ED";
  else if(ws[sz]=='y'&&ws[sz-1]=='l')uk=uk+"LY";
  else if(ws[sz]=='s')uk=uk+"S";
  else if(ws[sz]=='t'&&ws[sz-1]=='s'&&ws[sz-2]=='e')uk=uk+"EST";
  else if(ws[sz]=='r'&&ws[sz-1]=='e')uk=uk+"ER";
  else if(ws[sz]=='n'&&ws[sz-1]=='o'&&ws[sz-2]=='i')uk=uk+"ION";
  else if(ws[sz]=='e'&&ws[sz-1]=='l'&&ws[sz=2]=='b')uk=uk+"BLE";
  else if(ws[sz]=='t'&&ws[sz-1]=='n'&&ws[sz=2]=='e')uk=uk+"ENT";
  else if(ws[sz]=='y'&&ws[sz-1]=='r'&&ws[sz=2]=='o')uk=uk+"ORY";
  else if(ws[0]=='u'&&ws[1]=='n')uk="UN"+uk;
  else if(ws[0]=='e'&&ws[1]=='m')uk="EM"+uk;
  else if(ws[sz]=='l'&&ws[sz-1]=='a')uk=uk+"AL";
    for(;;i++){
      if(i==sz)break;
      if(ws[i]=='-'){
	uk=uk+"-";
	break;
      }
      else if(ws[i]=='.'){
	uk=uk+".";
	break;
      }
    }
  if(uk=="UNK")cerr<<ws<<endl;
  uk="*"+uk+"*";
  stIndex ans=Tree::tst->toIndex(uk);
  return ans;
}
*/
