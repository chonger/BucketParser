#include "Unkify.h"
#include "Tree.h"

bool
isUnk(stIndex w)
{
  string ws=Tree::tst->toString(w);
  int sz=ws.size();
  if(sz<5)return false;
  if(ws[0]!='*'||ws[sz-1]!='*')return false;
  return true;
}

stIndex
unkify(stIndex w)
{
  string ws=Tree::tst->toString(w);
  string uk="UNK";
  int sz=ws.size()-1;
  int i=0;
  if(isupper(ws[0]))uk="C"+uk;
  if(isdigit(ws[0])&&isdigit(ws[sz]))uk=uk+"N";
  else if(sz<=2){}
  else if(ws[sz]=='g'&&ws[sz-1]=='n'&&ws[sz-2]=='i')uk=uk+"ING";
  else if(ws[sz]=='d'&&ws[sz-1]=='e')uk=uk+"ED";
  else if(ws[sz]=='y'&&ws[sz-1]=='l')uk=uk+"LY";
  else if(ws[sz]=='s'&&ws[sz-1]=='u'&&ws[sz-2]=='o')uk=uk+"OUS";
  else if(ws[sz]=='s')uk=uk+"S";
  else if(ws[sz]=='t'&&ws[sz-1]=='s'&&ws[sz-2]=='e')uk=uk+"EST";
  else if(ws[sz]=='r'&&ws[sz-1]=='e')uk=uk+"ER";
  else if(ws[sz]=='n'&&ws[sz-1]=='o'&&ws[sz-2]=='i')uk=uk+"ION";
  else if(ws[sz]=='e'&&ws[sz-1]=='l'&&ws[sz=2]=='b')uk=uk+"BLE";
  else if(ws[sz]=='t'&&ws[sz-1]=='n'&&ws[sz=2]=='e')uk=uk+"ENT";
  else if(ws[sz]=='y'&&ws[sz-1]=='r'&&ws[sz=2]=='o')uk=uk+"ORY";
  else if(ws[0]=='u'&&ws[1]=='n')uk="UN"+uk;
  else if(ws[0]=='e'&&ws[1]=='m')uk="EM"+uk;
  else if(ws[sz]=='l'&&ws[sz-1]=='a')uk=uk+"AL";
    for(;;i++){
      if(i==sz)break;
      if(ws[i]=='-'){
	uk=uk+"-";
	break;
      }
      else if(ws[i]=='.'){
	uk=uk+".";
	break;
      }
    }

  uk="*"+uk+"*";
  stIndex ans=Tree::tst->toIndex(uk);
  return ans;
}

