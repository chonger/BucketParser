#ifndef PARSER_H
#define PARSER_H

#include "ECArgs.h"
#include "ECString.h"
#include "Tree.h"
#include "Rule.h"
#include "Unkify.h"
#include <assert.h>
#include "Trellis.h"
#include "CellClos.h"
#include <set>


#define OUTSIDE 1
#define INSIDE 0
#define VITERBI 2
#define NUMPARSERS 3

class Constit;
typedef vector<Constit*> ConstitV;

class CPair
{
 public:
  CPair() :l(NULL),r(NULL),p(NULL){}
  CPair(Constit* lh,Constit* rh,Rule* r) :l(lh),r(rh),p(r){}
  CPair(const CPair& cp) :l(cp.l),r(cp.r),p(cp.p){}
  Constit* l;
  Constit* r;
  Rule* p;
};

typedef vector<CPair> CPairs;
typedef CPairs::iterator CPairsIter;

class Constit
{
public:
  Constit(): label(0),inside(0),outside(0),mu(0),io(0),delta(0),rhs1(NULL),rhs2(NULL){}
    Constit(stIndex lab): label(lab),inside(0),outside(0),mu(0),io(0),delta(0),
    rhs1(NULL),rhs2(NULL){}
  Constit(const Constit& c):label(c.label),inside(c.inside),
    outside(c.outside),mu(c.mu),io(c.io),delta(c.delta),rhs1(c.rhs1),rhs2(c.rhs2),cpairs(c.cpairs){}
  stIndex label;
  double inside;
  double outside;
  double mu;
  double io;
  double delta;
  Constit* rhs1;
  Constit* rhs2;
  CPairs cpairs;
  Tree* makeTree(Indicies& unked,int& unkpos);
  void w(ostream& os){
    os<<Tree::ntst->toString(label);
    os<<" "<<inside<<" "<<outside<<" "<<inside*outside<<endl;
  }
};

typedef map<stIndex,Constit> Constits;
typedef Constits::iterator ConstitsIter;
typedef ConstitV::iterator ConstitVIter;
typedef set<Constit*> ConstitSet;

class Cell
{
public:
  Constit wrd;
  Constits constits;
};

class Parser;
typedef vector<Parser*> Parsers;

class ParserInfo
{
 public:
  ParserInfo(int pn,Rules* rule,string parf):pnum(pn),mode(0),rules(rule)
    {if(pn==NUMPARSERS-1)mode=VITERBI;
      if(pn>0)readPars(parf);
    }
  int pnum;
  int mode;
  Rules* rules;
  void readPars(string parf);
  stIndex pars[MAXNUMWORDS];
};

typedef vector<ParserInfo*> ParserInfos;
  
class Parser
{
public:
  Parser(Sent& snt, ParserInfo* pi,Parser* prevP,double ctf);
  Tree* parse(Sent& snt);
  Cell cells[MAXLEN+1][MAXLEN+1];
  void outside(int st,int ln, Cell& cell);
  void outsideUnary(Cell& cell);
  void outsideBinary(Cell& cell);
  void fill(int st,int ln, Cell& cell);
  void fillBinary(int st,int ln,int k,Cell& cell, int outsidep);
  void fillUnary(int st,int len,Cell& cell);
  void unkem(Indicies& wrds);
  void fixHeads();
  void countCs();
  void score(Tree* t,int& st);
  void earlyscore(Tree* t,int& st);
  bool ctfOk(int st,int nd,stIndex sti);
  void show();
  void prune();
  Rules* rules(){return pi->rules;}
  int& mode(){return pi->mode;}
  Tree* makeTree(Constit& c,Indicies& unked,int& unkpos,int& debg);
  stIndex par(stIndex lab){return pi->pars[lab];}
  static void init(ECString s);
  static void cleanup();
  static Parser* parse(vector<string>);
  double prob(string c, int st, int nd);
  double topInside;
  double ctf;
  int pnum;
  int len;
  Tree* result;
  Indicies words;
  Sent sent;
  static ECString dataDir;
  static long int numcs[2];
  static long int bucks[2][16];
  ParserInfo* pi;
  Parser* prevParser;
};

#endif
