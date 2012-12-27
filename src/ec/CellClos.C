#include "CellClos.h"
#include "Tree.h"
#include <assert.h>
#include "Trellis.h"
#include "Unkify.h"
#include "Parser.h"

#define igNumber 0
#define FUDGE 30

void
CellClos::
findTAnses(Tree* t,TAnses& tas,int& pos)
{
  if(!t->subtrees->subtrees){
    pos++;
    return;
  }
  int startpos=pos;
  for(Tree* p=t->subtrees;p;p=p->sibling)findTAnses(p,tas,pos);
  if(pos>(startpos+1)){
    tas[startpos].vals[STYPE]=true;
    tas[pos].vals[ETYPE]=true;
  }
}

void
CellClos::
fillIntPs(int d,IntPs& itps,int sp,Sent&sent)
{
  /* is sp-2 , so if sp=1 then we need stopWord for pos -1 */
  stIndex prv=0;
  stIndex prv2=0;
  stIndex prvw=0;
  stIndex prvw2=0;
  for(int i=(sp-(SFPOS/2)); i<(sp+(SFPOS/2));i++){
    stIndex sti=stopSym;
    stIndex stiw=stopWord;
    if(i>=0&&i<sent.sz){
      stiw=sent.unks[i];
      sti=sent.tags1[i];
    }
    Word& w=Trellis::words[stiw];
    int*  weiw=w.weights[d];
    int* aw=&(weiw[i-sp+(WFPOS/2)]);
    itps.push_back(aw);

    State& tag=Trellis::states[sti];
    int* wei= tag.weights[d];
    int* a1=&(wei[i-sp+(SFPOS/2)]);
    itps.push_back(a1);

    if(i>(sp-(SFPOS/2))){
      State& s2=Trellis::states2[prv][sti];
      int* wei2=s2.weights[d];
      int* a2=&(wei2[i-sp+(SFPOS/2)]);
      itps.push_back(a2);
      Wweights& ws=w.wweights[prvw];
      int* weiw=ws.w[d];
      int* aw2=&(weiw[i-sp+(SFPOS/2)]);
      itps.push_back(aw2);
      if(i>(sp+1-(SFPOS/2))){
	Wweights& s3=Trellis::states2[prv2][prv].w3s[sti];
	int* wei3=s3.w[d];
	int* a3=&(wei3[i-sp+(SFPOS/2)]);
	itps.push_back(a3);
	if(Trellis::words[prvw2].cnt>2000000){
	  Wweights2& ws3=ws.w3s[prvw2];
	  int* iw3=&(ws3.w[d][i-sp+(SFPOS/2)]);
	  itps.push_back(iw3);
	}
      }
    }
    prv2=prv;
    prv=sti;
    prvw2=prvw;
    prvw=stiw;
  }
}

  
void
CellClos::
procSpace(int sp,Sent& sent)
{
  for(int d=0;d<2;d++){
    IntPs weights;
    /* weights are the pointers to the weights supporting no/yes
       starting (if d=0) phrasal constit l>1 */
    int merits=0;

      fillIntPs(d,weights,sp,sent); 
      for(IntPsIter ipI=weights.begin();ipI!=weights.end();ipI++){
	merits+=**ipI;
      }

      //TAnses& anses=sent.anses;
      //TAns& tans=anses[sp];

      bool correct=false;//tans.vals[d];
    bool guessed=1;
    if(correct==0)gld++;
    bool ans=true;
    if(merits>0){
      guessed=0;
      if(merits>=FUDGE){
	tot++;
	ans=false;
	if(correct==0)totC++;
      }
    }
    else if(0>merits){
      guessed=1;
    }
    else if(doTraining){
      guessed =!correct;
    }
    //cerr<<"PS "<<sp<<" "<<ans<<" "<<merits<<" "<<sent.tags1[sp]<<" "<<sent.yield[sp]<<endl;
    if(d==0){
      TAns ta;
      sent.guessed.push_back(ta);
    }
    sent.guessed.back().vals[d]=ans;
    //cerr<<"GG "<<sp<<" "<<d<<" "<<ans<<endl;
    if(!doTraining)continue;
    if(guessed==correct)continue;

      IntPs& ws=weights;  
      int toadd= (!correct)? 1:-1;
      for(IntPsIter ipI=ws.begin();ipI!=ws.end();ipI++){
	(**ipI)+=toadd;
      }

  }
}

int
CellClos::
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

bool
hasnzf(int ws[][WFPOS])
{
  for(int nd=0;nd<NUMDSH;nd++){
    for(int l=0;l<WFPOS;l++){
      int v=ws[nd][l];
      if(v>0){
	return true;
      }
    }
  }
  return false;
}

void
CellClos::
writeFeats()
{
  int nd,l;
  for(int w=OFFSET+1;w<(Tree::tst->lword-1);w++){
    Word& word=Trellis::words[w];
    int numw2s=0;
    for(Wmap::iterator wmI=word.wweights.begin();
	wmI!=word.wweights.end();wmI++){
      stIndex w2=wmI->first;
      Wweights& ww=(wmI->second);
      if(!hasnzf(ww.w))continue;
      numw2s++;
    }
    cout<<Tree::tst->toString(w)<<" "<<word.cnt<<" "<<numw2s<<"\n\t";
    for(nd=0;nd<NUMDSH;nd++){
      for(l=0;l<WFPOS;l++)	cout<<word.weights[nd][l]<<" ";
    }
    cout<<endl;
    for(Wmap::iterator wmI=word.wweights.begin();
	wmI!=word.wweights.end();wmI++){
      stIndex w2=wmI->first;
      Wweights& ww=(wmI->second);
      if(!hasnzf(ww.w))continue;
      cout<<"\t"<<Tree::tst->toString(w2)<<" "<<ww.w3s.size()<<" ";
      for(int nd=0;nd<NUMDSH;nd++){
	for(int l=0;l<WFPOS;l++)  cout<<ww.w[nd][l]<<" ";
      }
      cout<<endl;
      for(Wmap2Iter wm2I=ww.w3s.begin();wm2I!=ww.w3s.end();wm2I++){
	Wweights2 ww2=wm2I->second;
	stIndex w3=wm2I->first;
	cout<<"\t\t"<<Tree::tst->toString(w3)<<" ";
	for(int nd=0;nd<NUMDSH;nd++){
	  for(int l=0;l<WFPOS;l++)  cout<<ww2.w[nd][l]<<" ";
	}
	cout<<endl;
      }
    }
  }
  for(int s=1;s<NUMSTATES;s++){
    State& st=Trellis::states[s];
    cout<<s<<"\t";
    for(nd=0;nd<NUMDSH;nd++)
      for(l=0;l<WFPOS;l++){
	cout <<st.weights[nd][l]<<" ";
      }
    cout<<endl;
  }
  for(int s=1;s<NUMSTATES;s++){
    for(int t=1;t<NUMSTATES;t++){
      int hasnz=0;
      State& st2=Trellis::states2[s][t];
      for(nd=0;nd<NUMDSH;nd++)	for(l=0;l<WFPOS;l++)
	if(st2.weights[nd][l]>0)hasnz=1;
      int numst2s=0;
      for(Wmap::iterator wmI=st2.w3s.begin();wmI!=st2.w3s.end();wmI++){
	Wweights& wws=wmI->second;
	stIndex st=wmI->first;
	if(!hasnzf(wws.w))continue;
	numst2s++;
      }
      cout<<s<<" "<<t<<" "<<numst2s<<"\t";
      for(nd=0;nd<NUMDSH;nd++){
	for(l=0;l<WFPOS;l++)	  cout<<st2.weights[nd][l]<<" ";
      }
      cout<<endl;
      for(Wmap::iterator wmI=st2.w3s.begin();wmI!=st2.w3s.end();wmI++){
	Wweights& wws=wmI->second;
	stIndex st=wmI->first;
	if(!hasnzf(wws.w))continue;

	cout<<"\t"<<st;
	for(nd=0;nd<NUMDSH;nd++){
	  for(l=0;l<WFPOS;l++){
	      cout<<" "<<wws.w[nd][l];
	  }
	}
	cout<<endl;
      }
      cout<<endl;
    }
  }
}

void
CellClos::
train(Sents& sents)
{
  for(iter=0; iter<MAXITER; iter++){
    cerr<<"ITER "<<iter<<endl;
    totC=tot=gld=0;
    for(int i=0;i<sents.size();i++){
      for(int sp=1;sp<(sents[i].sz);sp++){
	procSpace(sp,sents[i]);
      }
    }
  }
}

void
CellClos::
closer(Sent& sent)
{
  for(int sp=1;sp<(sent.sz);sp++){
    procSpace(sp,sent);
  }
  TAns ta;
  sent.guessed.push_back(ta);
  sent.guessed.back().vals[1]=1;
}


void
CellClos::
readFeats()
{
    ifstream ifs((Parser::dataDir +"weights.txt").c_str());
    assert(ifs);
  int nd,l;
  for(int w=OFFSET+1;w<(Tree::tst->lword-1);w++){
    string ckw;
    ifs>>ckw;
    //printf("!!!%s!\n",ckw.c_str());
    assert(ckw==Tree::tst->toString(w));
    int ct;
    ifs>>ct;
    int numsubs=0;
    ifs>>numsubs; 
    assert(numsubs>=0);
    Word& word=Trellis::words[w];
    word.cnt=ct;
    for(nd=0;nd<NUMDSH;nd++){
      for(l=0;l<WFPOS;l++){
	int w;
	ifs>>w;
	word.weights[nd][l]=w;

      }
    }
    for(int ns=0;ns<numsubs;ns++){
      string wrd;
      ifs>>wrd;
      stIndex w2=Tree::tst->toIndex(wrd);

      Wweights ww;
      int wwct;
      ifs>>wwct;
      for(int nd=0;nd<NUMDSH;nd++){
	for(int l=0;l<WFPOS;l++)  ifs>>ww.w[nd][l];
      }
      word.wweights[w2]=ww;
      for(int ctww=0;ctww<wwct;ctww++){
	string w3str;
	ifs>>w3str;
	stIndex w3idx=Tree::tst->toIndex(w3str);
	Wweights2 ww2;
	for(int nd=0;nd<NUMDSH;nd++){
	  for(int l=0;l<WFPOS;l++)  ifs>>ww2.w[nd][l];
	}
	ww.w3s[w3idx]=ww2;
      }
    }
  }
  for(int s=1;s<NUMSTATES;s++){
    State& st=Trellis::states[s];
    int cs;
    ifs>>cs;
    assert(cs==s);
    for(nd=0;nd<NUMDSH;nd++)
      for(l=0;l<WFPOS;l++){
	ifs>>st.weights[nd][l];
      }
  }
  for(;;){
    int s,t;
    int suscnt=0;
    ifs>>s;
    if(!ifs)break;
    assert(s>0);
    ifs>>t;
    assert(t>0);
    ifs>>suscnt;
    assert(suscnt>=0);
    State& st2=Trellis::states2[s][t];
    for(nd=0;nd<NUMDSH;nd++)
      for(l=0;l<WFPOS;l++){
	ifs>>st2.weights[nd][l];
      }
    for(int i=0;i<suscnt;i++){
      Wweights wws;
      int st;
      ifs>>st;
      assert(st>0);
      for(nd=0;nd<NUMDSH;nd++){
	for(l=0;l<WFPOS;l++){
	  ifs>>wws.w[nd][l];
	  }
	}
      st2.w3s[st]=wws;
    }
  }
}
