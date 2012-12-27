
#include "Parser.h"

int
main()
{
  Parser::init("/home/chonger/BucketParser/src/ec/DATA/");
  vector<string> vstr;
  vstr.push_back("I");
  vstr.push_back("ate");
  vstr.push_back("lunch");
  vstr.push_back("."); 
  Parser* p=Parser::parse(vstr);
  cout<<p->prob("NP",0,1)<<endl;
  cout<<p->prob("NP",1,2)<<endl;
  cout<<p->prob("S",0,4)<<endl;
  delete p;
  exit(0);
}
/* output I get
0.999917
0.000103939
0.994861
*/
