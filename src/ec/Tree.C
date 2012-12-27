#include "Tree.h"
#include <fstream>
#include <sys/resource.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

static pthread_mutex_t stlock = PTHREAD_MUTEX_INITIALIZER;
int Tree::toLower=0;

int Tree::singleLineReader=0; //set to 0 to read prettyprinted trees.
SymbolTable* Tree::ntst=NULL;
SymbolTable* Tree::tst=NULL;
Tree::Tree(Tree& other):
	label(other.label),
	head(other.head),
	subtrees(NULL),
	sibling(NULL),
	parent(NULL),
	pred(NULL),fnVals(other.fnVals)
{
	if(other.subtrees)
	{
		subtrees = new Tree(*other.subtrees);

		for(Tree* st = subtrees;
			st != NULL;
			st = st->sibling)
		{
			st->parent = this;
		}
	}

	if(other.sibling)
	{
		sibling = new Tree(*other.sibling);
	}
}

void
fix0(Tree* t)
{
  stIndex nonI=Tree::ntst->toIndex("-NONE-");
  if(t->label!=nonI)return;
  Tree* p=t->subtrees;
  if(!p)return;
  stIndex zeroI=Tree::tst->toIndex("0");
  if(p->label!=zeroI)return;
  p->label=Tree::tst->toIndex("*0*");
}

void
Tree::yield(vector<stIndex>& vsti)
{
  if(!subtrees)vsti.push_back(label);
  else{
    Tree* p;
    for(p=subtrees;p;p=p->sibling) p->yield(vsti);
  }
}

void
Tree::w(ostream& os){write(os);os<<endl;}

void
Tree::
writeYield(ostream& os)
{
  if(!subtrees)os<<tst->toString(label)<<" ";
  else{
    Tree* p;
    for(p=subtrees;p;p=p->sibling) p->writeYield(os);
  }
}

void
Tree::
write(ostream& os)
{
  if(!subtrees) os << tst->toString(label);
  else{
    Tree* p;
    os<< "(" << ntst->toString(label);
    for(p=subtrees;p;p=p->sibling){
      os<<" ";
      p->write(os);
    }
    os<<")";
  }
}

void
skipspaces(istream& is,int& cp,string& buf)
{
  int c;
  for(;;){
    if(cp>=buf.size()){
      if(Tree::singleLineReader)return;
      buf="";
      cp=0;
      getline(is,buf);
      if(!is)return;
    }
    c=buf[cp++];
    if(!isspace(c))break;
  }
  cp--;
}

void
readLabelStuff(Tree* t, char *s, int brkpt,int n)
{
  int k=0,i;
  char tmp[256];
  for(i=brkpt;i<n;i++){
    if(s[i]=='#'){
      int startnum;
      k=startnum=i+1;

      if(s[i+1]!='-'){
	for(k=i+1;k<n;k++){
	  if(s[k]=='-')break;
	  tmp[k-brkpt-1]=s[k];
	}
	tmp[k-brkpt-1]='\0';
	startnum=++k;
	//fprintf(stderr,"OB %s %i\n",tmp,n);
      }
      else{
	k++; 
	startnum++;
      }
      for(;k<n;k++){
	if(s[k]== '~')break;
	tmp[k-startnum]=s[k];
      }
      tmp[k-startnum]='\0';
    }
  }
}

void
readLabelOpen(Tree* t, char *s, int brkpt,int n)
{
  int k=0,i;
  char tmp[256];
  for(i=brkpt;i<n;i++){
    if(s[i]=='#'){
      int startnum;
      k=startnum=i+1;
      if(s[i+1]!='-'){
	for(k=i+1;k<n;k++){
	  if(s[k]=='-')break;
	  tmp[k-brkpt-1]=s[k];
	}
	tmp[k-brkpt-1]='\0';
	startnum=++k;
      }
      else{
	k++; 
	startnum++;
      }
    }
  }
}

stIndex
readlabel(istream& is,Tree* t,int& cp,string& buf)
{
  int	 c, brkpt=0,n=0, i=0;
  char s[256];
  for(c= buf[cp++] ;(cp<=buf.size() && !isspace(c) && c!='(' && c!=')');c=buf[cp++]){
    s[n++]=c;
    assert(n<256);
  }
  cp--;
  brkpt=n;
  for (i=1; i<n-1; i++){		// don't look at 1st or last char 
    if (strchr(CATSEP,s[i])) {		// if s[i] is in CATSEP
      brkpt= i;
      break;
    }
  }
  readLabelOpen(t,s,brkpt,n);
  s[brkpt] = '\0';
  pthread_mutex_lock(&stlock);
  stIndex sti=Tree::ntst->toIndex(s);
  pthread_mutex_unlock(&stlock);
  t->label=sti;
  return sti;
}

stIndex
readGWord(istream& is,int& cp, string& buf)
{
  char	 s[MAXLABELLEN];

  int	 n=0;
  char   c;
  for(c=buf[cp++];(cp<=buf.size() &&!isspace(c) && (c!='(') && (c!=')'));c=buf[cp++]){
    if(Tree::toLower)c=tolower(c);
    s[n++] = c;
    assert(n<MAXLABELLEN-1);		/* leave space for '\0' */
  }  
  s[n] = '\0';
  cp--;
  if(n == 0) return 0;
  pthread_mutex_lock(&stlock);
  stIndex ans = Tree::tst->toIndex((string)s);
  pthread_mutex_unlock(&stlock);
  return ans;
}

void
addParent(Tree* tr,Tree* par)
{
  if(!tr)return;
  if(par) tr->parent=par;
  Tree* p;
  for(p=tr->subtrees;p;p=p->sibling)addParent(p,tr);
}

Tree*
Tree::
make(istream& is)
{
  int cp=0;
  string buf;
  if(!singleLineReader)skipspaces(is,cp,buf);
  else{
    getline(is,buf);
  }
  if(!is) return NULL;
  if(buf.empty()&&singleLineReader){
    Tree* ns = new Tree(NEWSTORY);
    return ns;
  }
  Tree *ans=readB(is,cp,buf);
  addParent(ans,NULL);
  return ans;
}

Tree*
Tree::
readB(istream& is,int& cp, string& buf)
{
  int 	c;
  Tree  *t,*p;
  skipspaces(is,cp,buf);

  c=buf[cp];
  cp++;
  switch (c) {

  case ')': 
    return(NULL); break;	/* empty tree */

  case '(':			/* nonterminal */
    skipspaces(is,cp,buf);
    t = new Tree(0);
    readlabel(is,t,cp,buf);
    t->subtrees = p = readB(is,cp,buf);
    while (p) {
      Tree* stree=readB(is,cp,buf);
      p->sibling = stree;
      p = p->sibling;
    }
    assert(t->label);
    fix0(t);
    return(t); break;
  
  default:		
    cp--;
    skipspaces(is,cp,buf);
    t = new Tree(0);
    t->sibling = t->subtrees=NULL;
    t->label = readGWord(is,cp,buf);
    if(t->label<=0){
      cerr<<buf<<endl;
      assert(t->label>0);
    }
    return(t);
  }}

  
int
Tree::
wPos()
{
  int pos=0;
  Tree* p;
  for(p=this;p->parent;p=p->parent){}
  if(p->wPos(this, pos))return pos;
  else return -1;
}
    
int
Tree::
wPos(Tree* ht,int& pos)
{
  Tree* p=subtrees;
  if(!p){
    pos++;
    return 0;
  }
  else if(this==ht) return 1;
  for(;p;p=p->sibling) if(p->wPos(ht,pos))return 1;
  return  0;
}
