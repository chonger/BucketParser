//.959373 //.959434  .276/.028 seconds
//.962214   .665/.032
#include <assert.h>
#include <ctype.h>
#include "Trellis.h"
#include "Unkify.h"
#include "Parser.h"
#define MAXNUMWORDS 100000
#define MULTFAC 300

stIndex rootSym;
stIndex stopSym;
stIndex stopWord;

Word Trellis::words[MAXNUMWORDS];
State Trellis::states[NUMSTATES];
State Trellis::states2[NUMSTATES][NUMSTATES];
FT Trellis::fudge=0.96;
int tot1=0;
int tot2=0;

stIndex enstate(stIndex fst, stIndex sec) {return fst*NUMSTATES+sec;}
stIndex secstate(stIndex bth){ return bth%NUMSTATES; }
stIndex fststate(stIndex bth){ return bth/NUMSTATES; }

void
Sent::
unkem()
{
  int sz=yield.size();
  for(int pos=0;pos<sz;pos++){
    stIndex w=yield[pos];
    if(w>=Tree::tst->lword){
      unked.push_back(w);
      unks.push_back(unkify(w));
    }
    else unks.push_back(w);
  }
}

void
Trellis::
readSyms()
{
  ifstream ntfs((Parser::dataDir + "nTerms.txt").c_str());
  for(;;){
    string nt;
    ntfs>>nt;
    if(!ntfs)break;
    Tree::ntst->toIndex(nt);
  }
}

void
Trellis::
init(SymbolTable& ntst, SymbolTable& tst)
{
  readSyms();
  rootSym=ntst.toIndex(string("Root"));
  stopWord=tst.toIndex(string("sToP"));
  stopSym=ntst.toIndex(string("STOP"));
  Trellis::readData();
  tst.lword=tst.size()+OFFSET;
}

Trellis::
Trellis(Sent& snt):sz(snt.sz)
{
  for(int i=0;i<NUMSTATES2;i++)tstates[0][i].mu=0.0;
  tstates[0][enstate(stopSym,stopSym)].mu=1.0;
  yield=snt.unks;
  forward();
  backward();
  answer2(multiresult);
  snt.idv=multiresult;
  result=answer();
  snt.tags1=result;
}

typedef vector< pair<stIndex, FT> >::iterator vpsfIter;

void Trellis::
readData()
{ 
  ifstream ifs((Parser::dataDir + "ndata.txt").c_str()); 
  assert(ifs);
  int statecnts[NUMSTATES];
  for(int i=1;i<NUMSTATES;i++)statecnts[i]=0;
  for(int i=1;i<NUMSTATES;i++){
    for(int j=1;j<NUMSTATES;j++){
      State& st2=states2[i][j];
      for(int k=1;k<NUMSTATES;k++){
	string dum;
	ifs>>dum;
	ifs>>dum;
	ifs>>dum;
	int ct;
	ifs>>ct;
	statecnts[i]+=ct;
	st2.ntag+=ct;
	st2.nntag[k]+=ct;
      }
    }
    for(int j=1;j<NUMSTATES;j++){
      State& st2=states2[i][j];
      int stnt2=st2.ntag;
      FT v=(FT)stnt2/(FT)statecnts[i];
      st2.prb2=MULTFAC*v;
      for(int k=1;k<NUMSTATES;k++){
	int stnnt2=st2.nntag[k];
	if(stnnt2>0)assert(stnt2>0);
	if(stnt2==0)continue;
	FT v2=(FT)stnnt2 /(FT)stnt2;
	assert(v2<=1.0);
	st2.prob2[k]=MULTFAC*v2;
      }
    }
  }
  for(;;){
    string ws;
    ifs>>ws;
    stIndex x=Tree::tst->toIndex(ws);
    Word& w=words[x];
    for(;;){
      int pos;
      ifs>>pos;
      if(!ifs)break;
      if(pos==9999999)break;
      int ct;
      ifs>>ct;
      FT wp=(FT)ct/statecnts[pos];
      assert(wp<=1);
      pair<stIndex,FT> pr(pos,wp);
      w.wprobs.push_back(pr);
      //cerr<<ws<<" "<<x<<" "<<pos<<" "<<w.probs[pos]<<endl;
    }
    if(!ifs)break;
  }
}

/* I eat .    sz=3, find pointer at (3,stop) */
Indicies
Trellis::
answer()
{
  Indicies rans;
  stIndex prev=-1;
  FT bst=0;
  
  for(int s=1;s<NUMSTATES;s++){
    stIndex sm1=enstate(s,stopSym);
    FT m=tstates[sz][sm1].mu;
    if(m>bst){
      bst=m;
      prev=sm1;
    }
  }
  
  rans.push_back(fststate(prev));
  for(int loc=sz;loc>1;loc--){
    int nxt=tstates[loc][prev].bestprv;
    rans.push_back(fststate(nxt));
    prev=nxt;
  }
  Indicies ans;
  for(Indicies::reverse_iterator rI=rans.rbegin();rI!=rans.rend();rI++){
    ans.push_back(*rI);
  }
  return ans;
}


void
Trellis::
answer2(IndiciesV& ans)
{
  stIndex nxt=0;
  for(int loc=0;loc<sz;loc++){
    Indicies dummy;
    ans.push_back(dummy);
    FT tgs[NUMSTATES];
    for(int i=1;i<NUMSTATES;i++)tgs[i]=0;
    stIndex w=wrd(loc);
    Word& word=words[w];
    int wm1=wrd(loc-1);
    Word& word1=words[wm1];
    for(vpsfIter vpsfI=word.wprobs.begin();vpsfI!=word.wprobs.end();vpsfI++){
      int j=vpsfI->first;
      FT wprob=vpsfI->second;
      FT mx=0;
      int bst=-1;
      for(vpsfIter vpsfI3=word1.wprobs.begin();vpsfI3!=word1.wprobs.end();
          vpsfI3++){
        int i=vpsfI3->first;
        stIndex st=enstate(i,j);
        FT v=tstates[loc][st].alpha;
        tgs[j]+=v;
      }
    }
    tot1++;
    multimap<FT,stIndex> mm;
    for(int i=1;i<NUMSTATES;i++){
      pair<FT,stIndex> pr(tgs[i],i);
      mm.insert(pr);
    }
    multimap<FT,stIndex>::reverse_iterator rI=mm.rbegin();
    FT totsofar=0;
    for(;rI!=mm.rend();rI++){
      FT incr=rI->first;
      ans.back().push_back(rI->second);
      totsofar+= incr;
      tot2++;
      //cerr<<loc<<" "<<incr<<" "<<rI->second<<" "<<ans.back().size()<<endl;
      if(totsofar>=fudge)break;
    }
  }
  //cerr<<tot2<<" "<<tot1<<" "<<(FT)tot2/(FT)tot1<<endl;
}


void
Trellis::
forward()
{
  for(int loc=0;loc<=sz;loc++){
    stIndex w=wrd(loc);
    Word& word=words[w];
    int wm2=wrd(loc-2);
    Word& word2=words[wm2];
    int wm1=wrd(loc-1);
    Word& word1=words[wm1];
    for(vpsfIter vpsfI=word.wprobs.begin();vpsfI!=word.wprobs.end();vpsfI++){
      int j=vpsfI->first;
      FT wprob=vpsfI->second;
      FT mx=0;
      int bst=-1;
      for(vpsfIter vpsfI3=word1.wprobs.begin();vpsfI3!=word1.wprobs.end();
	  vpsfI3++){
	FT alph=0;
	int i=vpsfI3->first;
	for(vpsfIter vpsfI2=word2.wprobs.begin();vpsfI2!=word2.wprobs.end();
	    vpsfI2++){
	  int h=vpsfI2->first;
	  FT mu=0;
	  FT al=0;
	  if(loc==0){
	    if(i==stopSym){
	      mu=1;
	      al=1;
	    }
	  }
	  else{
	    mu=tstates[loc-1][enstate(h,i)].mu;
	    al=tstates[loc-1][enstate(h,i)].alpha;
	  }
	  FT ans =mu;
	  assert(j<50);
	  FT prb= states2[h][i].prob2[j];
	  if(prb==0)	prb=.000001;
	  ans*= prb;
	  alph+=prb*al;
	  if(ans>mx){
	    mx=ans;
	    bst=enstate(h,i);
	  }
	}
       tstates[loc][enstate(i,j)].alpha=(alph*wprob);
       //cerr<<"AL "<<loc<<" "<< i<<" "<<j<< " = " <<alph*wprob<<endl;
      }
      stIndex ns=enstate(secstate(bst),j);  
      mx*= wprob;
      tstates[loc][ns].mu=mx;
      tstates[loc][ns].bestprv=bst;
    }
  }
}
void
Trellis::
backward()
{
  FT sumAlpha=0;
  for(int i=1;i<NUMSTATES;i++) sumAlpha+=tstates[sz][enstate(i,stopSym)].alpha;
  sumAlpha*=MULTFAC/2;
  for(int loc=sz;loc>=0;loc--){
    stIndex w=wrd(loc-1);
    Word& word=words[w];
    int wm2=wrd(loc+1);
    Word& word2=words[wm2];
    int wm1=wrd(loc);
    Word& word1=words[wm1];
    for(vpsfIter vpsfI=word.wprobs.begin();vpsfI!=word.wprobs.end();vpsfI++){
      int j=vpsfI->first;
      for(vpsfIter vpsfI3=word1.wprobs.begin();vpsfI3!=word1.wprobs.end();
          vpsfI3++){
        FT beta=0;
        int i=vpsfI3->first;
        for(vpsfIter vpsfI2=word2.wprobs.begin();vpsfI2!=word2.wprobs.end();
            vpsfI2++){
          int h=vpsfI2->first;
          FT probofw2gh=vpsfI2->second;
          FT bet=0;
          if(loc==sz){
            if(i==stopSym){
              bet=1.0/sumAlpha;
            }
          }
          else{
            bet=tstates[loc+1][enstate(i,h)].beta;
          }
	  assert(h<50);
          FT prb= states2[j][i].prob2[h];
          if(prb==0)    prb=.000001;
          beta+=prb*bet*probofw2gh;
        }
        tstates[loc][enstate(j,i)].beta=beta;
        FT alph=tstates[loc][enstate(j,i)].alpha;
        tstates[loc][enstate(j,i)].alpha*=beta;
        //cerr<<"BB "<<loc<<" "<<j<<" "<<i<<" "<<beta<<" "<< alph*beta<<endl;
      }
    }
  }
}
