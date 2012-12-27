#ifndef CELLCLOS_H
#define CELLCLOS_H

#include "SymbolTable.h"
#include "Trellis.h"

typedef vector<int*> IntPs;
typedef IntPs::iterator IntPsIter;

#define MAXITER 30

class CellClos
{
public:
  CellClos():doTraining(0),totC(0),tot(0),gld(0),iter(0){}
  void  findTAnses(Tree* t,TAnses& tas,int& pos);
  void fillIntPs(int d,IntPs& itps,int sp,Sent&sent);
  void procSpace(int sp,Sent& sent);
  int score(Indicies& gold, Indicies& ans);
  bool sufficient(int ws[][WFPOS]);
  static void readFeats();
  void writeFeats();
  void setFeats();
  void train(Sents& sents);
  void closer(Sent& sent);
  int doTraining;
  int totC;
  int tot;
  int gld;
  int iter;
  vector<bool> results[2];
};

#endif
