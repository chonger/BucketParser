
#ifndef TRELLIS_H
#define TRELLIS_H

#define FT double
#include "ECArgs.h"
#include "ECString.h"
#include "Tree.h"

#define NUMSTATES 49
#define MAXLEN 100
#define MAXNUMWORDS 100000
#define NUMSTATES2 NUMSTATES*NUMSTATES

#define NUMDSH 2
#define STYPE 0
#define ETYPE 1
#define SFPOS 4
#define WFPOS 4

extern stIndex rootSym;
extern stIndex stopSym;
extern stIndex stopWord;

class Wweights2
{
 public:
  Wweights2()
    {for(int i=0;i<NUMDSH;i++)for(int j=0;j<WFPOS;j++)w[i][j]=0;}
  Wweights2(const Wweights2& ww)
    {for(int i=0;i<NUMDSH;i++)for(int j=0;j<WFPOS;j++)w[i][j]=ww.w[i][j];}
  int w[NUMDSH][WFPOS];
};

typedef map<stIndex, Wweights2> Wmap2;
typedef Wmap2::iterator Wmap2Iter;

class Wweights
{
 public:
  Wweights()
    {for(int i=0;i<NUMDSH;i++)for(int j=0;j<WFPOS;j++)w[i][j]=0;}
  Wweights(const Wweights& ww)
    {for(int i=0;i<NUMDSH;i++)for(int j=0;j<WFPOS;j++)w[i][j]=ww.w[i][j];}
  int w[NUMDSH][WFPOS];
  Wmap2 w3s;
};


typedef map<stIndex, Wweights> Wmap;
typedef Wmap::iterator WmapIter;

class Word
{
public:
  Word():cnt(0){
    for(int i=0;i<NUMDSH;i++)for(int j=0;j<WFPOS;j++)weights[i][j]=0;
  }
  int cnt;
  vector< pair<stIndex,FT> > wprobs;
  Wmap wweights;
  Wweights* fww(stIndex i){
    WmapIter wmI=wweights.find(i);
    if(wmI!=wweights.end())return &(wmI->second);
    return NULL;
  }
  int weights[NUMDSH][WFPOS];
};
      


class State
{
public:
  State():cnt(0),prb2(0),ntag(0)
  {
    for(int i=0;i<50;i++){
      nntag[i]=0;
    }
    for(int i=0;i<NUMDSH;i++)for(int j=0;j<SFPOS;j++)weights[i][j]=0;
  }
  int cnt;
  FT prb2;
  int ntag;
  int nntag[50];
  FT prob2[50];
  int weights[NUMDSH][SFPOS];
  Wmap w3s;
};

class TAns
{
public:
  TAns(){for(int i=0;i<NUMDSH;i++)vals[i]=false;}
  bool vals[NUMDSH];
};

typedef vector<TAns> TAnses;
typedef vector<TAnses> TAnsesV;

class Sent
{
public:
  Sent(Indicies& y):sz(y.size()),yield(y)
    {TAns ta;
      guessed.push_back(ta);
      guessed.back().vals[0]=1;
      unkem();}
  Sent(Indicies& y,TAnses& tas):sz(y.size()),yield(y),anses(tas)
    {TAns ta;
      guessed.push_back(ta);
      guessed.back().vals[0]=1;
      unkem();}
  Sent(const Sent& s): sz(s.sz),yield(s.yield),tags1(s.tags1),
    unked(s.unked),unks(s.unks),idv(s.idv),anses(s.anses),guessed(s.guessed){}
  int sz;
  Indicies yield;
  Indicies tags1;
  Indicies unked;
  Indicies unks;
  IndiciesV idv;
  vector<bool> open;
  TAnses anses;
  TAnses guessed;
  int pos(int pos,int wh){
    return pos+wh;
  }
  stIndex tags(int i){
    if(i<0)return stopSym;
    if(i>=sz)return stopSym;
    else return tags1[i];
  }
  void wordPos(int i,Indicies& idis){
    if(i<0)idis.push_back(stopSym);
    else idis=idv[i];
  }
  void unkem();
};  

typedef vector<Sent> Sents;

class TState
{
public:
  TState():bestprv(0),mu(0),alpha(0),beta(0){}
  stIndex bestprv;
  FT mu;
  FT alpha;
  FT beta;
};

class Trellis
{
public:
  Trellis(Sent& s);
  void forward();
  void backward();
  void compForward(int i, stIndex j,FT wprb);
  Indicies answer();
  void answer2(IndiciesV& iv);
  stIndex wrd(int loc){
    if(loc<0)return stopWord;
    else if(loc>=sz)return stopWord;
    else return yield[loc];
  }
  static void init(SymbolTable& ntst,   SymbolTable& tst);
  static void readData();
  static FT fudge;
  int sz;
  Indicies yield;
  Indicies result;
  IndiciesV multiresult;
  TState tstates[MAXLEN+2][NUMSTATES2];
  static Word words[MAXNUMWORDS];
  static State states[NUMSTATES];
  static State states2[NUMSTATES][NUMSTATES];
  static void readSyms();
};

#endif
