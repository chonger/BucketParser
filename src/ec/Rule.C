
#include "Rule.h"
#include <sstream>

bool
binarized(stIndex label)
{
  string labstring=Tree::ntst->toString(label);
  for(int i=0;i<labstring.size();i++){
    if(labstring[i]=='_')return true;
  }
  return false;
}

string
headed(string labstring)
{
  string ans;
  bool ish=false;
  for(int i=0;i<labstring.size();i++){
    if(labstring[i]=='@'){
      ish=true;
      break;
    }
    else ans+=labstring[i];
  }
  if(!ish)return "";
  else return ans;
}


void
Rule::
w(ostream& os)
{
  os<<prob<<" "<<Tree::ntst->toString(lhs)<<" ";
  if(rhs1>Tree::tst->offset){
    os<<"-> ";
    os<<Tree::tst->toString(rhs1);
  }
  else{
    os<<"--> ";
    os<<Tree::ntst->toString(rhs1);
    if(rhs2)os<<" "<<Tree::ntst->toString(rhs2)<<" ";
  }
  os<<endl;
}

Rule*
Rules::
findB(Rule* r)
{
  RuleMap1& rm1=brules[r->rhs1];
  RuleMap2& rm2=rm1[r->rhs2];
  Rule* rr =rm2[r->lhs];
  if(!rr){
    Rule* nr= new Rule(*r);
    rm2[r->lhs]=nr;
    return nr;
  }
  return rr;
}
  
Rule*
Rules::findU(Rule* r)
{
  RuleMap1& rm1=urules[r->typ];
  RuleMap2& rm2=rm1[r->rhs1];
  Rule* rr =rm2[r->lhs];
  if(!rr){
    Rule* nr= new Rule(*r);
    rm2[r->lhs]=nr;
    urules[2][r->lhs][r->rhs1]=nr;
    return nr;
  }
  return rr;
}
  
Rule*
Rules::
find(Rule* r)
{
  if(r->typ<2)return findU(r);
  else return findB(r);
}

void
Rules::
readr(ifstream& ifs)
{
  assert(ifs);
  int sz=Tree::tst->offset;
  int denom[OFFSET];
  for(int i=0;i<OFFSET;i++)denom[i]=0;
  int nrules=0;
  for(;;){
    stringstream ss;
    string ln;
    getline(ifs,ln);
    ss<<ln;
    if(!ifs)break;
    int cnt;
    ss>>cnt;
    string arrow,lhs, rhs1,rhs2;
    stIndex lhsi, rhs1i, rhs2i;
    ss>>lhs;
    lhsi=Tree::ntst->toIndex(lhs);
    ss>>arrow;
    ss>>rhs1;
    rhs2i=0;
    int typ=0;
    short bt=0;
    if(arrow=="->"){
      rhs1i=Tree::tst->toIndex(rhs1);
      string hh=headed(lhs);
      if(hh!=""){
	stIndex orgh=Tree::ntst->toIndex(hh);
	//indicate that the word index rhs1i, when seen with orgh(IN@of) should
	//also be considered lhsi (IN).
	fixheadMap[rhs1i][orgh]=lhsi;
      }
    }
    else{
      rhs1i=Tree::ntst->toIndex(rhs1);
      typ=1;
      if(!ss.eof()){
	typ=2;
	ss>>rhs2;
	rhs2i=Tree::ntst->toIndex(rhs2);
	if(binarized(lhsi))bt=ZEDLHS;
      }
    }
    Rule r(typ,lhsi,rhs1i,rhs2i);
    if(rhs2i==0&& rhs1i==lhsi){
      continue;
    }
    r.lhstyp=bt;
    r.cnt=cnt;
    denom[lhsi]+=cnt;

    Rule* rl=find(&r);
    nrules++;
  }
  for(int i=0;i<OFFSET;i++){
    RuleMap1& rules1=brules[i];
    for(RuleMap1::iterator rmI1=rules1.begin();rmI1!=rules1.end();rmI1++){
	RuleMap2& rules2=rmI1->second;
	for(RuleMap2::iterator rmI2=rules2.begin();rmI2!=rules2.end();rmI2++){
	  Rule* r=rmI2->second;
	  assert(r);
	  r->prob = (double)(r->cnt)/(double)denom[r->lhs];
	}
    }
  }
  for(int i=0;i<2;i++){
    RuleMap1& rules1=urules[i];
    for(RuleMap1::iterator rmI1=rules1.begin();rmI1!=rules1.end();rmI1++){
      RuleMap2& rules2=rmI1->second;
      for(RuleMap2::iterator rmI2=rules2.begin();rmI2!=rules2.end();rmI2++){
	Rule* r=rmI2->second;
	assert(r);
	r->prob = (double)(r->cnt)/(double)denom[r->lhs];
      }
    }
  }
}
      
    
	     
