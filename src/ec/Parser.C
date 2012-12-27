#include "Parser.h"
#include <math.h>
#include "Trellis.h"

#define REG 0
#define ZED 1
#define CENTERC '/'
long int Parser::numcs[2]={0,0};
long int Parser::bucks[2][16]={{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};
int debug = 0;
void
wi(stIndex i){cerr<<" "<<Tree::ntst->toString(i);}

double
Parser::
prob(string c,int st, int nd)
{
  assert(nd<MAXLEN);
  Cell& cl=cells[st][nd];
  stIndex lab=Tree::ntst->toIndex(c);
  ConstitsIter ci= cl.constits.find(lab);
  if(ci==cl.constits.end())return 0;
  return ci->second.io;
}

stIndex
parentIndex(stIndex lab)
{
  if(!lab)return 0;
  string labs=Tree::ntst->toString(lab);
  string nstring;
  for(int i=0;i<labs.size();i++){
    assert(labs[i]!='|');
    if(labs[i]=='&'||
       labs[i]=='+'||
       labs[i]=='^'||
       (i>0&&labs[i]=='@'))break;
    else nstring+=labs[i];
  }
  stIndex dep=Tree::ntst->toIndex(nstring);
  return dep;
}

int
isCenter(stIndex lab)
{
  string labs=    Tree::ntst->toString(lab);
  int sz=labs.size();
  if(sz<3)return false;
  for(int i=1;i<labs.size();i++)if(labs[i]==CENTERC)return 1;
  return 0;
}

void
deParent(Tree* t)
{
  if(!t->subtrees)return;
  for(Tree* p=t->subtrees;p;p=p->sibling)deParent(p);
  string lab=Tree::ntst->toString(t->label);
  string nstring;
  for(int i=0;i<lab.size();i++){
    assert(lab[i]!='|');
    if(lab[i]=='&'||
       lab[i]=='+'||
       lab[i]=='^'||
       (i>0&&lab[i]=='@'))break;
    else nstring+=lab[i];
  }
  stIndex dep=Tree::ntst->toIndex(nstring);
  t->label=dep;
}

Parser::
Parser(Sent& snt,ParserInfo* pari,Parser* prevp,double ctf):
  sent(snt),pi(pari), prevParser(prevp),ctf(ctf)
{
  mode()=INSIDE;
  fixHeads();
  words=snt.unks;
  result=NULL;
  len=words.size();
  //because we never use ln==0, we need to make the array size MAXLEN+1
  for(int ln=1;ln<=len;ln++){
    for(int start=ln-1;start>=0;start--){
      assert(start>=0);
      assert(start<len);
      assert(len<=MAXLEN);
      fill(start,ln,cells[start][ln]);
    }
  }
  Constit& top=cells[0][len].constits[rootSym];
  if(mode() ==VITERBI){
    if(top.mu==0)return;
    int up=0;
    int debg=0;
    result=makeTree(top,sent.unked,up,debg);
    deParent(result);
    assert(result);
    return;
  }

  topInside=top.inside;
  if(topInside==0)return;
  top.outside=1.0/top.inside;
  //cerr<<"TO="<<top.inside<<endl;
  for(int ln=len;ln>=1;ln--){
    for(int st=0;st<ln;st++){
      outside(st,ln,cells[st][ln]);
    }
  }
  if(pi->pnum==3)countCs();

}

void
Parser::
fillUnary(int st,int nd,Cell& cell)
{
  RuleMap1& rm1=rules()->urules[0];
  if(nd-st==1){
    stIndex wsti=cell.wrd.label;
    Indicies ids;
    sent.wordPos(st,ids);//sets ids;
    bool seen=false;
    for( IndiciesIter iI=ids.begin();iI!=ids.end();iI++){
      stIndex lhs=*iI;
      assert(lhs>0);
      Rule* r=rm1[wsti][lhs];
      if(!r)continue;
      if(r->lhs==rootSym&&(st!=0||nd!=(len-1)))continue;
      Constit con(lhs);
      seen=true;
      con.inside=r->prob;
      con.rhs1=&cell.wrd;
      con.mu=r->prob;
      cell.constits[lhs]=con;
      assert(con.label!=0);
      //cerr<<"FU "<<st<<" ";     con.w(cerr);
    }
    assert(seen);
  }
  RuleMap1& rm1u=rules()->urules[1];
  ConstitV curAdd;
  for(ConstitsIter cI=cell.constits.begin();cI!=cell.constits.end();cI++){
    Constit& con=cI->second;
    con.delta=con.inside;
    curAdd.push_back(&con);
  }
  ConstitSet alreadyDone;
  while(!curAdd.empty()){
    ConstitV toAdd;
    for(ConstitVIter vcI=curAdd.begin();vcI!=curAdd.end();vcI++){
      Constit* cadd=*vcI;
      double delta=cadd->delta;
      cadd->delta=0;
      RuleMap2& rm2=rm1u[cadd->label];
      for(RuleMap2::iterator rm2I=rm2.begin();rm2I!=rm2.end();rm2I++){
	Rule* r=rm2I->second;
	assert(r);
	assert(r->lhs>0);

	if(prevParser)	  if(!ctfOk(st,nd,r->lhs))	    continue;

	Constit& oC=cell.constits[r->lhs];
	if(mode()==VITERBI){
	  double mu=r->prob*cadd->mu;
	  if(oC.mu>0)	    assert(oC.label!=0);
	  if(oC.mu<mu){
	    oC.mu=mu;
	    oC.label=r->lhs;
	    oC.rhs1=cadd;
	    oC.rhs2=NULL;
	    toAdd.push_back(&oC);
	    assert(oC.label!=0);
	  }
	}
	else{
	  //oC is parent, cadd is child;
	  double inside=r->prob*delta;
	  assert(mode()==INSIDE);
	  if(oC.inside<(5*inside))	    toAdd.push_back(&oC);
	  oC.inside+=inside;
	  oC.delta+=inside;
	  oC.label=r->lhs;
	  //so this cp says we built oC from cadd using rule r.
	  if(alreadyDone.find(cadd)==alreadyDone.end()){
	    CPair cp(cadd,NULL,r);
	    oC.cpairs.push_back(cp);
	  }
	  //cerr<<oC.cpairs.size()<<" "; 	  oC.w(cerr);
	  //cerr<<"\t";r->w(cerr);
	}
	assert(oC.label);
      }
      if(mode()!=VITERBI)alreadyDone.insert(cadd);
    }
    curAdd=toAdd;
  }
}

void
Parser::
fill(int st,int nd, Cell& cell)
{
  if(nd-st==1){
    Constit con(words[st]);
    con.inside=1;
    con.mu=1;
    cell.wrd=con;
    //cerr<<"WD "<<Tree::tst->toString(con.label)<<endl;
  }
  if(st==0||sent.guessed[st].vals[0]==1){
    for(int k=st+1;k<nd;k++)
      if(k==(nd-1)||sent.guessed[nd].vals[1]==1)
	fillBinary(st,nd,k,cell,mode()==VITERBI?VITERBI:INSIDE);
  }
  fillUnary(st,nd,cell);
}

void
Parser::
outside(int st,int nd, Cell& cell)
{
  //cerr<<"TT "<<st<<" "<<nd<<" "<<cell.constits.size()<<endl;
  outsideUnary(cell);
  outsideBinary(cell);
  for(ConstitsIter cI=cell.constits.begin();cI!=cell.constits.end();cI++){
    Constit& ci=cI->second;
    ci.io=ci.inside*ci.outside;
    if(ci.io>2){cerr<<"BAD "<<st<<" "<<nd<<" ";ci.w(cerr);}
  }
}

bool 
Parser::
ctfOk(int st,int nd,stIndex lab)
{
  Cell& pc=prevParser->cells[st][nd];
  stIndex pidx=par(lab);
  if(!pidx){
    cerr<<st<<" "<<nd<<" "<<Tree::ntst->toString(lab)<<endl;
    assert(pidx);
  }
  ConstitsIter cI=pc.constits.find(pidx);
  if(cI==pc.constits.end()){
    return false;
  }
  Constit& con=cI->second;
  if(debug)cerr<<"\t"<<"CTF "<<con.inside*con.outside<<endl;
  if(debug&&lab<OFFSET)cerr<<"\t"<<Tree::ntst->toString(lab)<<endl;
  if(con.io<ctf)return false;
  return true;
}

void
Parser::
fillBinary(int st,int end,int k,Cell& cell, int modefb)
{
  Cell& rhs1=cells[st][k];
  Cell& rhs2=cells[k][end];
  for(ConstitsIter cI=rhs1.constits.begin();
      cI!=rhs1.constits.end();cI++){
    Constit& r1c=cI->second;
    RuleMap1& rm1=rules()->brules[r1c.label];
    for(RuleMap1::iterator rm1I=rm1.begin(); rm1I!=rm1.end();rm1I++){
      RuleMap2& rm2=rm1I->second;
      stIndex rhs2i=rm1I->first;
      ConstitsIter c2I=rhs2.constits.find(rhs2i);
      if(c2I==rhs2.constits.end())continue;
      Constit& r2c=c2I->second;
      for(RuleMap2::iterator rm2I=rm2.begin();rm2I!=rm2.end();rm2I++){
	Rule* r=rm2I->second;
	stIndex lhs=r->lhs;
	if(r->lhs==rootSym&&st==0&&end==(len-1))continue;
	if(prevParser&&!ctfOk(st,end,lhs))continue;
	double rp=r->prob;
	Constit& con=cell.constits[lhs];
	if(modefb==INSIDE){
	  double m1=r1c.inside;
	  double m2=r2c.inside;
	  double inside=rp*m1*m2;
	  con.label=lhs;
	  con.inside+=inside;
	  CPair cp(&r1c,&r2c,r);
	  con.cpairs.push_back(cp);
	}
	else{
	  con.label=lhs;
	  double m1=r1c.mu;
	  double m2=r2c.mu;
	  double mu=rp*m1*m2;
	  if(con.mu<mu){
	    con.mu=mu;
	    con.rhs1=&r1c;
	    con.rhs2=&r2c;
	  }
	  assert(con.label!=0);
	  assert(mu>0);
	  //if(st==5&&end==5&&pi->pnum==4)
	  if(false)
	    cerr<<st<<" "<<k<<" "<<end<<" "<<Tree::ntst->toString(lhs)<<" "<<con.mu<<" "<<m1<<" "<<m2<<" "<<rp<<" "<<con.cpairs.size(
)<<endl; 
	}
      }
    }
  }
}

void
Parser::
fixHeads()
{
  IndiciesVIter ivI=sent.idv.begin();
  //go through all the words of the sentence;
  for(IndiciesIter iI=sent.yield.begin();iI!=sent.yield.end();iI++){
    StiMapMapIter stimm=rules()->fixheadMap.find(*iI);
    //if the word is an annotated one
    if(stimm!=rules()->fixheadMap.end()){
      Indicies& iv=(*ivI);
      StiMap& stim=stimm->second;
      //go through allowed postags
      for(IndiciesIter iI=iv.begin();iI!=iv.end();iI++){
	stIndex org=*iI;
	stIndex anpos=stim[org];
	//this clobbers e.g., IN, with IN@of ???
	if(anpos)(*iI)=anpos;
      }
    }
    ivI++;
  }
}



void
Parser::
countCs()
{
  for(int ll=1;ll<=len;ll++){
    for(int st=ll-1;st>=0;st--){
      Cell& cl=cells[st][ll];
      for(ConstitsIter cI=cl.constits.begin();
	  cI!=cl.constits.end();cI++){
	Constit& cc=cI->second;
	if(cc.inside>0)numcs[0]++;
	double io=cc.inside*cc.outside;
	assert(io>=0);
	if(io<ctf)bucks[0][15]++;
	if(io==0)bucks[0][14]++;
	else{
	  int b=-log2(io);
	  if(b>14)b=14;
	  bucks[0][b]++;
	}
      }
    }
  }
}

stIndex
child(Tree* t,int wh)
{
  stIndex ans=t->label;
  if(!t->subtrees||!t->subtrees->subtrees)return ans;
  if(wh==0){
    if(t->label==rootSym)return rootSym;
    return Tree::ntst->toIndex("S");
  }
  else if(wh==1){
    string splhs=Tree::ntst->toString(t->label);
    if(splhs[0]=='N')ans=Tree::ntst->toIndex("NP");
    else if(splhs=="VP")ans=Tree::ntst->toIndex("VP");
    else if(splhs=="PP")ans=Tree::ntst->toIndex("PP");
    else if(t->label==rootSym)ans=rootSym;
    else ans=Tree::ntst->toIndex("S");
  }
  return ans;
}

void
Parser::
show()
{
  for(int l=len;l>0;l--)
    for(int s=0;s<l;s++)
      for(ConstitsIter cI=cells[s][l].constits.begin()
	    ;cI!=cells[s][l].constits.end();cI++){
	Constit& c=cI->second;
	if(c.io>.1&&!binarized(c.label)){
	  cerr<<s<<" "<<l<<" ";c.w(cerr);
	}
      }
}

void
Parser::
score(Tree* t,int& st)
{
  if(pi->pnum<3)return earlyscore(t,st);
  int s=st;
  if(!t->subtrees->subtrees){
    st++;
    return;
  }
  Tree* p=t->subtrees;
  for(;p;p=p->sibling)    score(p,st);
  numcs[1]++;
  int nd=st;
  Constits& ccs=cells[s][nd].constits;
  double maxio=0;
  for(ConstitsIter cI=ccs.begin();cI!=ccs.end();cI++){
    Constit& con=cI->second;
    //lab is the complex label in chart, par is it's simplified version
    stIndex lab=con.label;
    if(!lab){
      cerr<<"ZERO LABEL "<<st<<" "<<nd<<" ";
      t->w(cerr);
    }
    stIndex par=parentIndex(lab);

    if(par==t->label){
      /*      
      cerr<<"SEE "<<s<<" "<<nd<<" "
	<<Tree::ntst->toString(par)<<" "
	<<Tree::ntst->toString(lab)<<" "<<Tree::ntst->toString(t->label)
	<<" "<<con.io<<endl;
      */
      double io=con.io;
      if(io>maxio)maxio=io;
    }
  }
  if(maxio<ctf)bucks[1][15]++;
  if(maxio==0)bucks[1][14]++;
  else{
    int b=-log2(maxio);
    if(b>14)b=14;
    bucks[1][b]++;
  }
}

void
Parser::
earlyscore(Tree* t,int& st)
{
  int s=st;
  if(!t->subtrees->subtrees){
    st++;
    return;
  }
  Tree* p=t->subtrees;
  for(;p;p=p->sibling)    earlyscore(p,st);
  int nd=st;
  stIndex sym=child(t,pi->pnum);
  if(!sym){
    t->w(cerr);
    assert(sym);
  }
  Constit& cc=cells[s][nd].constits[sym];
  numcs[1]++;
  double io=cc.inside*cc.outside;
  //cerr<<s<<" "<<nd<<" "<<io<<" "<<cc.label<<" "<<sym<<" "; t->w(cerr);
  if(io<ctf)bucks[1][15]++;
  if(io==0)bucks[1][14]++;
  else{
    int b=-log2(io);
    if(b>14)b=14;
    bucks[1][b]++;
  }
}


Tree*
Constit::
makeTree(Indicies& unked,int& unkpos)
{
  if(!rhs1){
    if(isUnk(label)){
      label=unked[unkpos++];
    }
    Tree* ans=new Tree(label);
    //ans->w(cerr);
    return ans;
  }
  Tree* t1=rhs1->makeTree(unked,unkpos);
  Tree* t2=NULL;
  if(rhs2)    t2=rhs2->makeTree(unked,unkpos);
  else if(isCenter(label)){
    return t1;
  }

  Tree* s=t1;
  for(;s->sibling;s=s->sibling){}
  s->sibling=t2;

  if(binarized(label)){
    return t1;
  }
  Tree* t3=new Tree(label);
  t3->subtrees=t1;
  return t3;
}

Tree*
Parser::
makeTree(Constit& c,Indicies& unked,int& unkpos,int& debg)
{
  int st=debg;
  stIndex label = c.label;
  if(!c.rhs1){
    if(isUnk(label)){
      label=unked[unkpos++];
    }
    Tree* ans=new Tree(label);
    debg++;
    //cerr<<"WW ";ans->w(cerr);
    return ans;
  }
  Tree* t1=makeTree(*c.rhs1,unked,unkpos,debg);
  Tree* t2=NULL;
  if(c.rhs2)    t2=makeTree(*c.rhs2,unked,unkpos,debg);
  else if(isCenter(label)){
    //cerr<<ctfOk(st,debg,label)<<" "<<Tree::ntst->toString(label)<<" "<<st<<" "<<debg<<" ";    t1->w(cerr);
    return t1;
  }

  Tree* s=t1;
  for(;s->sibling;s=s->sibling){}
  s->sibling=t2;

  if(binarized(label)){
    //cerr<<ctfOk(st,debg,label)<<" "<<Tree::ntst->toString(label)<<endl;
    return t1;
  }
  Tree* t3=new Tree(label);
  t3->subtrees=t1;
  //cerr<<ctfOk(st,debg,label)<<": "<<Tree::ntst->toString(label)<<" "<<st<<" "<<debg<<" ";t3->w(cerr);
  return t3;
}

void
ParserInfo::
readPars(string pf)
{
  ifstream ifs(pf.c_str());
  if(!ifs){
    cerr<<pf<<endl;
    assert(ifs);
  }
  for(;;){
    string chi;
    string par;
    ifs>>chi;
    if(!ifs)break;
    ifs>>par;
    stIndex ch=Tree::ntst->toIndex(chi);
    assert(ifs);
    stIndex pa=Tree::ntst->toIndex(par);
    pars[ch]=pa;
  }
}

void
Parser::
outsideBinary(Cell& cell)
{
  for(ConstitsIter cI=cell.constits.begin();cI!=cell.constits.end();cI++){
    Constit& con=cI->second;
    if(debug)con.w(cerr);
    if(con.outside==0)continue;
    for(CPairsIter cpI=con.cpairs.begin();cpI!=con.cpairs.end();cpI++){
      CPair& cp=*cpI;
      if(!cp.r)continue;
      double ir=cp.r->inside;
      double il=cp.l->inside;
      double bi=con.outside*cp.p->prob;
      cp.r->outside+=bi*il;
      cp.l->outside+=bi*ir;
      //cerr<<"OB "; con.w(cerr);cerr<<"\t";      cp.l->w(cerr);       cerr<<"\t";cp.r->w(cerr);
    }
  }
}

void
Parser::
outsideUnary(Cell& cell)
{
  ConstitV toAdd;
  for(ConstitsIter cI=cell.constits.begin();cI!=cell.constits.end();cI++){
    Constit& con=cI->second;
    if(con.outside==0)continue;
    con.delta=con.outside;
    toAdd.push_back(&con);
  }
  while(!toAdd.empty()){
    ConstitV newAdd;
    for(ConstitVIter cvI =toAdd.begin();cvI!=toAdd.end();cvI++){
      Constit& parc=**cvI;
      double delta=parc.delta;
      parc.delta=0;
      for(CPairsIter cpI=parc.cpairs.begin();cpI!=parc.cpairs.end();cpI++){
	CPair& cp=*cpI;
	if(cp.r)continue;
	Constit& chic=*(cp.l);
	double newo=delta*cp.p->prob;
	double oldo=chic.outside;
	chic.outside+=newo;
	chic.delta+=newo;
	if(oldo<5*newo){
	  newAdd.push_back(&chic);
	  //cerr<<"OU ";chic.w(cerr);
	  //cerr<<"\t";parc.w(cerr);
	}
	assert(chic.outside*chic.inside<1.8);
      }
    }
    //cerr<<"NI"<<endl;
    toAdd=newAdd;
    newAdd.empty();
  }
}
      
      
void
Parser::
prune()
{
  for(int ln=len;ln>=1;ln--){
    for(int st=0;st<ln;st++){
      Indicies toRemove;
      Cell& cl=cells[st][ln];
      for(ConstitsIter cI=cl.constits.begin();cI!=cl.constits.end();cI++){
	Constit& con=cI->second;
	if(con.io<ctf)	  toRemove.push_back(con.label);
	else{
	  CPairs copyInto;
	  for(CPairsIter cpI=con.cpairs.begin();cpI!=con.cpairs.end();cpI++){
	    CPair& cp=*cpI;
	    if(cp.l->label>OFFSET){
	      copyInto.push_back(cp);
	      continue;
	    }
	    if(cp.l->io<ctf)continue;
	    else if(cp.r&&(cp.r->io<ctf))continue;
	    copyInto.push_back(cp);
	  }
	  /*
	  if(copyInto.empty()){
	    assert(con.io>=ctf);
	    cerr<<"EE "<<st<<" "<<ln<<" ";
	    con.w(cerr);
	  }
	  */
	  con.cpairs=copyInto;
	}
      }
      for(IndiciesIter iI=toRemove.begin();iI!=toRemove.end();iI++){
	stIndex ci=*iI;
	cl.constits.erase(ci);
      }
    }
  }
}
