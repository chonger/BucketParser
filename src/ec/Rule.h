
#ifndef RULE_H
#define RULE_H
#include "ECString.h"
#include "SymbolTable.h"
#include "Tree.h"

class Rule;

typedef map< stIndex, Rule*> RuleMap2;
typedef  map< stIndex, RuleMap2> RuleMap1;
typedef  map< stIndex, RuleMap1> RuleMap ;
typedef map<stIndex, stIndex> StiMap;
typedef map<stIndex, StiMap> StiMapMap;
typedef StiMapMap::iterator StiMapMapIter;

#define WORD 0
#define UNARY 1
#define BINARY 2
#define ZEDLHS 1

class Rule
{
public:
  Rule(Tree* t):cnt(0),prob(0),lhs(t->label),rhs1(t->subtrees->label){
    Tree* st2=(t->subtrees->sibling);
    if(st2) rhs2=st2->label;
    else rhs2=0;
  }
  void w(ostream& os);
  Rule(const Rule& r):typ(r.typ),lhstyp(r.lhstyp),lhs(r.lhs),rhs1(r.rhs1),
    rhs2(r.rhs2),cnt(r.cnt),prob(r.prob){}
  Rule(int otyp,stIndex olhs,stIndex orhs,stIndex orhs2):prob(0),typ(otyp),
      lhs(olhs),rhs1(orhs),rhs2(orhs2){}
  short typ;
  short lhstyp;
  int cnt;
  double prob;
  stIndex lhs;
  stIndex rhs1;
  stIndex rhs2;
};


class Rules
{
 public:
  Rules(string ifs):parent(0),child(0){
    ifstream ifstr(ifs.c_str());
    readr(ifstr);}
    ~Rules(){
        for(int i=0;i<OFFSET;++i) {
            for(RuleMap1::iterator j=brules[i].begin();j != brules[i].end();++j) {
                RuleMap2& rm2 = j->second;
                for(RuleMap2::iterator k=rm2.begin();k != rm2.end();++k) {
                    delete k->second;
                }
            }
        }
        for(int i=0;i<2;++i) {
            for(RuleMap1::iterator j=urules[i].begin();j != urules[i].end();++j) {
                RuleMap2& rm2 = j->second;
                for(RuleMap2::iterator k=rm2.begin();k != rm2.end();++k) {
                    delete k->second;
                }
            }
        }
    }
  int parent;
  int child;
  RuleMap1 brules[OFFSET];
  RuleMap1 urules[3];
  StiMapMap fixheadMap;
  Rule* findU(Rule* r);
  Rule* find(Rule* r);
  Rule* findB(Rule* r);
  void readr(ifstream& ifs);
};


bool binarized(stIndex label);
void readRules(ifstream& ifs);

#endif
