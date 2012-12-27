
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "Parser.h"
#include <sstream>
#include "Trellis.h"
#include "Tree.h"

#define numThreads 4

ifstream* ifs=NULL;
int ntrees=0;
ofstream* ofs=NULL;
ofstream* pofs=NULL;
#define CTF 0.0001
static pthread_mutex_t readlock = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t writelock = PTHREAD_MUTEX_INITIALIZER;
int totB=0;
int totB2=0;
ParserInfos parserinfos;

Parser*
Parser::
parse(vector<string> swrds)
{
  if(swrds.size()>MAXLEN)    return NULL;

  Indicies wrds;
  for(vector<string>::iterator sI=swrds.begin();sI!=swrds.end();sI++)
    wrds.push_back(Tree::tst->toIndex(*sI));
  Sent sent(wrds);
  Trellis trel(sent);
  CellClos clClos;
  clClos.closer(sent);
  Parser* prev=NULL;
  Parsers parsers;
  for(int i=0;i<NUMPARSERS;i++){
    Parser* pa = new Parser(sent,parserinfos[i],prev,CTF);
    parsers.push_back(pa);
    prev=pa;
  }
  for(int i=i;i<NUMPARSERS-1;i++)delete parsers[i];
  return prev;
}

void
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

string
itoa(int i)
{
  stringstream ss;
  ss<<i;
  return ss.str();
}

void
readParsers()
{
  //  Rules* rls = new Rules(Parser::dataDir + "benrules.txt");
  //parserinfos.push_back(new ParserInfo(0,rls,""));
  
  string rulepth=Parser::dataDir + "rules";
  string prpath=Parser::dataDir + "pars";
  for(int i=0;i<NUMPARSERS;i++){
    string rulepath=rulepth+itoa(i)+".txt";
    Rules* rls = new Rules(rulepath);
    string parspath="";
    if(i>0)parspath=prpath+itoa(i-1)+itoa(i)+".txt";
    parserinfos.push_back(new ParserInfo(i,rls,parspath));
  }
  
}

ECString Parser::dataDir;

void
Parser::
init(ECString s)
{
  Parser::dataDir = s;
  SymbolTable* ntst = new SymbolTable;
  SymbolTable* tst = new SymbolTable;
  tst->offset=OFFSET;
  Tree::ntst=ntst;
  Tree::tst=tst;
  Trellis::init(*ntst,*tst);
  rootSym=ntst->toIndex("Root");
  CellClos::readFeats();
  readParsers();
}

void Parser::cleanup() {
  for(size_t i=0;i<parserinfos.size();++i) {
    ParserInfo* pi = parserinfos[i];
    delete pi->rules;
    delete pi;
  }
}
  
