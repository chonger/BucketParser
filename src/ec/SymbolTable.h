
#ifndef STBL_H
#define STBL_H

#include "assert.h"
#include <map>
#include <vector>
#include "ECString.h"
#include <iostream>

typedef unsigned int stIndex;
typedef unsigned int usi;
typedef map<string,usi> SymMap;
typedef SymMap::iterator SymMapIter;

typedef vector<stIndex> Indicies;
typedef vector<Indicies> IndiciesV;
typedef Indicies::iterator IndiciesIter;
typedef IndiciesV::iterator IndiciesVIter;

#define OFFSET 70000

class SymbolTable
{
 public:
  SymbolTable(int offset1=0):offset(offset1),fword(0),lword(0){sarray.push_back("");}
  string   toString(stIndex i) {
    i-=offset;
    if(!(i>0&&i<sarray.size())){
      cerr<<"Bad index "<<i<<endl;
      assert(0);
    }
    else return sarray[i];
  }
  int size(){ return sarray.size(); }
  stIndex toIndex(string s){
    SymMapIter smi=smap.find(s);
    if(smi== smap.end()){
      int nn=sarray.size();
      sarray.push_back(s);
      smap[s]=nn;
      return nn+offset;
    }
    else return smi->second+offset;
  }
  stIndex isThere(string s){
    SymMapIter smi=smap.find(s);
    if(smi== smap.end())return 0;
    else return smi->second+offset;
  }
  /*
  stIndex toindex(string s){
    SymMapIter smi=smap.find(s);
    if(smi== smap.end()) return 0;
    else return smi->second+offset;
  }
  */
  int offset;
  stIndex fword;
  stIndex lword;
  SymMap smap;
  vector<string> sarray;
};

#endif
