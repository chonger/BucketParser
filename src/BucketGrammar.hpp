#ifndef BUCKETGRAMMAR
#define BUCKETGRAMMAR 1

#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <google/sparsetable>
#include <google/dense_hash_map>
#include <google/dense_hash_set>
#include <set>
#include <boost/functional/hash.hpp>
#include "Chart.hpp"
using namespace std;
using boost::hash;

typedef google::dense_hash_map<unsigned int,set<unsigned int>,hash<unsigned int>,equal_to<unsigned int> > Umap;
typedef google::dense_hash_map<unsigned int,vector<pair<unsigned int,unsigned int> >,hash<unsigned int>,equal_to<unsigned int> > Lmap;
typedef google::dense_hash_map<pair<unsigned int, unsigned int>,set<unsigned int>,hash<pair<unsigned int , unsigned int> >,equal_to<pair<unsigned int, unsigned int> > > Bmap;
typedef google::dense_hash_map<unsigned int,double > I2Dmap;
typedef google::dense_hash_map<unsigned int,unsigned int > I2Imap;
typedef google::dense_hash_map<string,unsigned int,hash<string>,equal_to<string> > S2Imap;
typedef google::dense_hash_map<string,set<unsigned int>,hash<string>,equal_to<string> > S2Setmap;
typedef google::dense_hash_map<vector<unsigned int>,unsigned int,hash<vector <unsigned int> >, equal_to<vector <unsigned int> > > V2Imap;



struct TreeNode {

    TreeNode(unsigned int sym_,vector<TreeNode*> k_) : sym(sym_), k(k_) {}
    
    ~TreeNode() {
        for(vector<TreeNode*>::iterator iter = k.begin();iter != k.end();++iter) {
            delete *iter;
        }
    }

    size_t size() {
        size_t r = 1;
        for(size_t i=0;i<k.size();++i) {
            r += k[i]->size();
        }
        return r;
    }

    vector<TreeNode*> nodes() {
        vector<TreeNode*> r;
        r.push_back(this);
        for(vector<TreeNode*>::iterator iter = k.begin();iter != k.end();++iter) {
            vector<TreeNode*> kBelow = (*iter)->nodes();
            for(vector<TreeNode*>::iterator jter = kBelow.begin();jter != kBelow.end();++jter) {
                r.push_back(*jter);
            }
        }
        return r;
    }

    unsigned int sym;
    vector<TreeNode*> k;
    
};


//makes the statement that a symbol (sym) has a given prob of appearing over span i,j
struct EvalItem {

    std::string sym;
    unsigned int i,j;
    double prob;
    
    EvalItem(std::string sym_, unsigned int i_, unsigned int j_, double prob_) : sym(sym_), i(i_), j(j_), prob(prob_) {
        
    }

    bool operator==(EvalItem& o) {
        return sym == o.sym && i == o.i && j == o.j;
    }
        
    
};

typedef google::dense_hash_map<pair<unsigned int, unsigned int>,vector<EvalItem*>,hash<pair<unsigned int , unsigned int> >,
                               equal_to<pair<unsigned int, unsigned int> > > EvalMap;

struct ParseTree {
    /**
    ParseTree(string s, S2Imap sym2base, unsigned int glueS) {
        stringstream ss(s);
        root = readNode(ss,sym2base,glueS);
        terms = getTerms(s);
    }
    */

    static int id;

    //int myID;
    
    ParseTree(string x) : root(NULL), terms(), hashv(x) {
        //myID = id++;
        //printf("!CREATE ID : %d\n",myID);
    }
    
    ParseTree(TreeNode* root_, vector<string> terms_) : root(root_), terms(terms_), hashv("") {
        //myID = id++;
        //        printf("CREATE ID : %d\n",myID);
    }

    void setHashString(string* syms) {
        vector<string>::iterator titer = terms.begin();
        hashv = getHString(root, titer, syms);
        //printf("TREE %s\n",hashv.c_str());
    }
    
    string getHString(TreeNode* n, vector<string>::iterator& titer, string* syms) {
        string ret = "(" + syms[n->sym];
        if(n->k.size() == 0) {
            ret += " " + *titer + ")";
            titer++;
            return ret;
        }
        for(size_t i=0;i<n->k.size();++i) {
            ret += " " + getHString(n->k[i],titer,syms);
        }
        ret += ")";
        return ret;
    }
    
    ~ParseTree() {
        //printf("DELETE ID : %d\n",myID);
        delete root;
    }

    size_t size() {
        return root->size();
    }

    size_t depth() {
        int ind = 0;
        return recDepth(root,ind);
    }
    
    size_t recDepth(TreeNode*& n, int& ind) {
        size_t ret = 0;
        if(n->k.size() == 0) {
            if(terms[ind] != "<>") {
                ret = 1;
            }
            ind++;
        } else {
            for(vector<TreeNode*>::iterator iter = n->k.begin();iter != n->k.end();++iter) {
                ret = max(ret,recDepth(*iter,ind)+1);
            }
        }
        return ret;
    }
    
    bool isPCFG() {
        return depth() == 1;
    }    
    
    vector<pair<TreeNode*,int> > nodes() {
        int i = 0;

        vector<pair<TreeNode*,int> > r;

        recNodes(root,i,r);
        
        return r;
    }

    vector<ParseTree*> rules() {
        vector<ParseTree*> ret;
        vector<pair<TreeNode*,int> > ns = nodes();
        for(size_t i=0;i<ns.size();++i) {
            vector<TreeNode*> ruleK;
            TreeNode* n = ns[i].first;
            int tI = ns[i].second;
            vector<string> t;
            if(n->k.size() == 0) {
                t.push_back(terms[tI]);
            } else {
                for(size_t j=0;j<n->k.size();++j) {
                    TreeNode* nn = n->k[j];
                    vector<TreeNode*> kk;
                    ruleK.push_back(new TreeNode(nn->sym,kk));
                    t.push_back("<>");
                }
            }
            ret.push_back(new ParseTree(new TreeNode(n->sym,ruleK),t));
        }
        return ret;
    }
    
    void recNodes(TreeNode* n, int& tInd, vector<pair<TreeNode*,int> >& ret) {

        if(n->k.size() == 0) {//terminal
            pair<TreeNode*,int> me =  make_pair(n,tInd++);
            ret.push_back(me);
        } else {
            for(size_t i=0;i<n->k.size();++i) {
                TreeNode* kn = n->k[i];
                recNodes(kn,tInd, ret);
            }
            pair<TreeNode*,int> me =  make_pair(n,-1);
            ret.push_back(me);
        }
    }
    
    //gets the gold standard EvalItems from a tree (when called on the root
    //returns the eval items and the number of terminals covered
    //TODO : syms is just the grammar symbols - could either change EvalItem to contain ints instead of strings, or make this a BucketGrammar member.
    std::pair<std::vector<EvalItem>,unsigned int> getEItems(TreeNode* node, string* syms) {
        
        std::vector<EvalItem> ret;
        
        if(node->k.size() == 0) {
            EvalItem e(syms[node->sym],0,0,1.0);
            ret.push_back(e);
            return make_pair(ret,1);
        } else {
            unsigned int tot = 0;
            for(size_t i=0;i<node->k.size();++i) {
                TreeNode* kk = node->k[i];
                pair< vector< EvalItem >, unsigned int> res = getEItems(kk,syms);
                for(size_t j=0;j<res.first.size();++j) {
                    EvalItem& ee = res.first[j];
                    ee.j += tot;
                    ret.push_back(ee);
                }
                tot += res.second;
            }
            
            EvalItem e(syms[node->sym],tot-1,0,1.0);
            ret.push_back(e);
            return make_pair(ret,tot);
        }     
    }
    
    TreeNode* root;
    vector<string> terms;
    string hashv;
};



struct ParseTreeEq {
  bool operator()(ParseTree* const& a, ParseTree* const& b ) const
  {
      return a->hashv == b->hashv;
  }
};


struct ParseTreeHash{
    size_t operator()(ParseTree* const& k) const{
        hash<string> hf;
        return hf(k->hashv);
    }
};

class BucketGrammar {
public:

    BucketGrammar() : adjoinP(NULL), syms(NULL) {

    }
    
    BucketGrammar(std::ifstream& ifs) : adjoinP(NULL), syms(NULL) {
        loadGrammar(ifs);
    }

    void loadGrammar(std::ifstream&ifs) {

        //HASMAP SETUP
        //these hashmaps are actual data members
        tsgP.set_empty_key(-1);
        prior.set_empty_key(-1);
        tagP.set_empty_key(-1);
        aClassMap.set_empty_key(-1);
        
        umap.set_empty_key(-1);
        umap.set_deleted_key(-2);

        leftlook.set_empty_key(-1);
        leftlook.set_deleted_key(-2);
        
        bmap.set_empty_key(make_pair(-2,-2));
        bmap.set_deleted_key(make_pair(-1,-1));

        preterms.set_empty_key("EMPTYKEY");
        preterms.set_deleted_key("DELETEDKEY");

        sym2base.set_empty_key("EMPTYKEY");
        sym2base.set_deleted_key("DELETEDKEY");
        
        //these hashmaps are only used in the transform (loading only)

        //this one maps terminals to the base symbols
        S2Setmap pt2base;
        pt2base.set_empty_key("EMPTYKEY");
        pt2base.set_deleted_key("DELETEDKEY");

        V2Imap goodman;
        vector<unsigned int> eKey,dKey;
        eKey.push_back(0);
        dKey.push_back(1);
        goodman.set_empty_key(eKey);
        goodman.set_deleted_key(dKey);

        //Start reading the grammar file here
        
        string line;
        getline(ifs,line);
        nSym = strtol(line.c_str(),NULL,0);
        //        printf("LINE - %s\n",line.c_str());
        //printf("%u symbols\n",nSym);

        syms = new string[2*nSym+1]; //save an extra spot for our special GLUE symbol
        //adjoinP = new double[nSym];

        for(int i=0;i<nSym;++i) {
            getline(ifs,line);
            string sym = line;
            //getline(ifs,line);
            //double prob = atof(line.c_str());
            //printf("%s - %f\n",sym.c_str(),prob);
            if(sym == "ROOT")
                rootSym = i;
            syms[i] = sym;
            syms[i+nSym] = sym+"*";
            sym2base[sym] = i;
            sym2base[sym + "*"] = i + nSym;
            //adjoinP[i] = prob;
        }
        syms[nSym*2] = "GLUE";

        glueS = nSym*2;

        
        goodmanIndex = nSym*2+1; //skip both base and foot nodes, and glue sym
        for(unsigned int i=0;i<goodmanIndex;++i) { //map the base/foot/glue nodes to themselves in the basemap
            canL.push_back(false);
            canR.push_back(false);
            baseSym.push_back(i);
        }

        //read TSG rules
        
        getline(ifs,line);
        nTSG = strtol(line.c_str(),NULL,0);
        //printf("%u TSG rules\n",nTSG);

        if(nTSG == 0)
            printf("WARNING : this line indicates no TSG rules - %s\n",line.c_str());
        
        for(int i=0;i<nTSG;++i) {
            if(i>0 && i % 100000 == 0)
                printf("processed %d so far\n",i);
            getline(ifs,line);
            if(allow(line)) {
                string s = line;
                unsigned int ruleRoot = addRule(line,sym2base,pt2base,goodman);
                getline(ifs,line);
                double prob = atof(line.c_str());
                if(prob > 1.0 || prob < 0.0) {
                    printf("%s -> %f : Illegal probability\n",s.c_str(),prob);
                    throw "FAIL";
                }

                tsgP[ruleRoot] = prob;
            } else {
                getline(ifs,line);
            }
        }

        //read TAG rules
        
        getline(ifs,line);
        nTAG = strtol(line.c_str(),NULL,0);
        //printf("%u TAG rules\n",nTAG);
        //nTAG = 0;

        
        for(int i=0;i<nTAG;++i) {
            if(i>0 && i % 100000 == 0)
                printf("processed %d so far\n",i);
            getline(ifs,line);
            if(allow(line)) {
                unsigned int ruleRoot = addRule(line,sym2base,pt2base,goodman);

                /**
                //DIAGNOSTIC!!!
                ParseTree pt(line,sym2base,nSym*2);
                string ty = pt.checkfoot(pt.root,nSym);                
                wrapping.insert(make_pair(ruleRoot,ty + "\t" + line));
                //DIAGNOSTIC!!!
                */
                
                getline(ifs,line);
                double prob = atof(line.c_str());
                tagP[ruleRoot] = prob;
            } else {
                getline(ifs,line);
            }
        }
        
        //printf("Transformed Grammar has %u nodes\n",goodmanIndex);

        getline(ifs,line);
        nAC = strtol(line.c_str(),NULL,0);

        adjoinP = new double[nAC];

        aKeys.set_empty_key("");
        
        for(int i=0;i<nAC;++i) {
            getline(ifs,line);
            aKeys.insert(make_pair(line,i));
            getline(ifs,line);
            double prob = atof(line.c_str());
            adjoinP[i] = prob;
        }
        
        getG();
        for(unsigned int i =nSym*2+1;i<goodmanIndex;++i) {
            unsigned int bS = baseSym[i];
            if(ptsyms.count(bS) == 0 && bS != nSym*2 && tagP.count(i) == 0) { //get rid of preterminals and GLUE symbol
                string k = getAKey(i);
                if(aKeys.count(k) > 0)
                    aClassMap.insert(make_pair(i,aKeys[k]));
            }
        }        
        
        ifs.close();
                
    }

    //SYMBOL AND LEFTMOST CHILD
    string getAKey(unsigned int s) {
        ParseTree* t = getTree(s);
        string os = toString(t);
        delete t;
        return os;
        //        return syms[baseSym[s]];// + " " + syms[baseSym[gKids[s][0]]];
    }
    
    virtual bool allow(string& s) {
        return true;
    }
    
    ~BucketGrammar() {
        if(adjoinP != NULL)
            delete[] adjoinP;
        adjoinP = NULL;
        if(syms != NULL)
            delete[] syms;
        syms = NULL;
        
    }

        //from a flat string representation of a ParseTree, get a vector of terminals
    //TODO : can be made a static utility func
    vector<string> getTerms(string& s) {
        vector<string> terms;
        stringstream ss(s);
        string w;
        while(ss.good()) {
            ss >> w;
            if(w[0] != '(') {
                string tok = w.substr(0,w.find_first_of(')'));
                if(tok != "<>") {
                    terms.push_back(tok);
                }
            }
        }
        return terms;
    }
    
    TreeNode* readNode(stringstream& ss) {
        char c;
        ss >> c;
        assert(c == '(');

        string sym = "";
        ss >> sym;
        unsigned int symI = sym2base[sym];
        
        vector<TreeNode*> kids;
        
        while(true) {
            c = ss.peek();
            if(c == ' ') { //it's a space
                ss.ignore(1); //ignore that space
            } else if(c == '(') { //nonterminal
                //get child node index
                TreeNode* kid = readNode(ss);
                kids.push_back(kid);                
            } else if (c == ')') {
                ss.ignore(1); //burn closing paren
                while(kids.size() > 2) { //try to add glue rule
                    TreeNode* r = kids.back();
                    kids.pop_back();
                    TreeNode* l = kids.back();
                    kids.pop_back();
                    unsigned int glueI = glueS;
                    vector<TreeNode*> nKids;
                    nKids.push_back(l);
                    nKids.push_back(r);
                    
                    kids.push_back(new TreeNode(glueI,nKids));
                }

                return new TreeNode(symI,kids);
                
            } else { //terminal
                string term = "";
                while(ss.peek() != ')') {
                    term += ss.get();
                }
                ss.ignore(1); //burn closing paren
                vector<TreeNode*> noKids;
                return new TreeNode(symI,noKids);
            }
        }
        return NULL;
    }

    unsigned int getRule(stringstream& ss, S2Setmap& pt2base, V2Imap& goodman) {
        char c;
        ss >> c;
        assert(c == '(');

        string sym = "";
        ss >> sym;
        unsigned int symI = sym2base[sym];
        //        printf("got symbol %s (%u)\n",sym.c_str(),symI);
        
        vector<unsigned int> kSyms;
        kSyms.push_back(symI);
        
        while(true) {
            c = ss.peek();
            if(c == ' ') { //it's a space
                ss.ignore(1); //ignore that space
            } else if(c == '(') { //nonterminal
                //get child node index
                unsigned int index = getRule(ss,pt2base,goodman);
                kSyms.push_back(index);                
            } else if (c == ')') {
                ss.ignore(1); //burn closing paren
                //      printf("Finished Rule with sym %s and %lu kids\n",sym.c_str(),kSyms.size()-1);
                unsigned int index = 0; 
                
                while(kSyms.size() > 3) { //try to add glue rule
                    //  printf("K = %lu\n",kSyms.size());
                    unsigned int r = kSyms.back();
                    kSyms.pop_back();
                    unsigned int l = kSyms.back();
                    kSyms.pop_back();

                    unsigned int glueI = nSym*2;
                    
                    vector<unsigned int> gluerule;
                    gluerule.push_back(glueI);//glue symbol
                    gluerule.push_back(l);
                    gluerule.push_back(r);

                    unsigned int glueindex = 0; 
                    V2Imap::iterator fter = goodman.find(gluerule);
                    if(fter == goodman.end()) { //new node
                        goodman[gluerule] = goodmanIndex;
                        baseSym.push_back(glueI);
                        glueindex = goodmanIndex;

                        canL.push_back(false);
                        canR.push_back(false);
                        
                        canL[l] = true;
                        canR[r] = true;

                        leftlook[l].push_back(make_pair(r,goodmanIndex));
                        
                        bmap[make_pair(l,r)].insert(goodmanIndex);
                        
                        goodmanIndex++;
                    } else { //seen it
                        glueindex = fter->second;
                    }
                    kSyms.push_back(glueindex);
                }

                V2Imap::iterator fter = goodman.find(kSyms);
                if(fter == goodman.end()) { //new node
                    goodman[kSyms] = goodmanIndex;
                    baseSym.push_back(symI);
                    index = goodmanIndex;

                    canL.push_back(false);
                    canR.push_back(false);
                    
                    goodmanIndex++;

                    if(kSyms.size() == 3) {
                        //add index -> l r
                        //    printf("BR : %u -> %u %u\n",index,kSyms[1],kSyms[2]);
                        canL[kSyms[1]] = true;
                        canR[kSyms[2]] = true;
                        leftlook[kSyms[1]].push_back(make_pair(kSyms[2],index));
                        bmap[make_pair(kSyms[1],kSyms[2])].insert(index);
                    } else { //one child
                        //add index -> k
                        umap[kSyms[1]].insert(index);
                        //printf("UR : %u -> %u\n",index,kSyms[1]);
                    }
                    
                } else { //seen it
                    index = fter->second;
                }
                                
                return index;
            } else { //terminal
                string term = "";
                while(ss.peek() != ')') {
                    term += ss.get();
                }
                //                printf("got terminal %s\n",term.c_str());
                ss.ignore(1); //burn closing paren

                unsigned int index = 0;
                if(term == "<>") { //nonterminal leaf
                    index = sym2base[sym];
                } else { //preterminal node
                    //get index and add
                    //printf("T:%s\n",term.c_str());
                    S2Setmap::iterator fter = pt2base.find(term);

                    ptsyms.insert(symI);
                    
                    if(fter != pt2base.end()) {//seen this terminal before
                        set<unsigned int>& bsymset = fter->second; //the base syms that have been seen to parse this terminal
                        
                        if(bsymset.find(symI) == bsymset.end()) { //never found this pterm rule before
                            preterms[term].insert(goodmanIndex);
                            index = goodmanIndex;
                            baseSym.push_back(symI);
                            ++goodmanIndex;
                            canL.push_back(false);
                            canR.push_back(false);
                        
                            bsymset.insert(symI);
                        } else {
                            //this is a small set...maybe an unideal implementation tho
                            set<unsigned int>& symset = preterms[term];
                            for(set<unsigned int>::iterator iter = symset.begin();iter != symset.end();++iter) {
                                if(baseSym[*iter] == symI)
                                    index = *iter;
                            }
                        }
                    } else {
                        preterms[term].insert(goodmanIndex);
                        index = goodmanIndex;
                        baseSym.push_back(symI);
                        canL.push_back(false);
                        canR.push_back(false);
                        ++goodmanIndex;
                 
                        set<unsigned int> bsymset;
                        bsymset.insert(symI);
                        pt2base[term] = bsymset;
                    }
                }
                //                printf("returning %u - %u\n",index,goodmanIndex);
                return index;
            }
        }
        
        return 1;
    }

    
    unsigned int addRule(string& s, S2Imap& sym2base, S2Setmap& pt2base, V2Imap& goodman) {
        stringstream ss;
        ss << s;
        return getRule(ss,pt2base,goodman);
    }

    I2Dmap tsgP;
    I2Dmap tagP;
    I2Dmap prior;
    vector<unsigned int> baseSym;
    S2Setmap preterms; //maps terminal strings to sets of goodman transformed preterms
    Umap umap;
    Bmap bmap;
    S2Imap aKeys; //adjunction keys
    S2Imap sym2base;
    
    vector<bool> canL;
    vector<bool> canR;
    Lmap leftlook;
    
    unsigned int nSym;
    unsigned int glueS;
    int nTSG;
    int nTAG;
    int nAC;
    I2Imap aClassMap;
    double* adjoinP;
    string* syms;
    unsigned int goodmanIndex;
    unsigned int rootSym;
    vector<vector<unsigned int> > gKids;
    vector<string> tKids;
    set<unsigned int> ptsyms;

    void makePrior() {
        getG();
        for(I2Dmap::iterator iter = tsgP.begin();iter != tsgP.end();++iter) {
            unsigned int rSym = iter->first;
            ParseTree* tree = getTree(rSym);
            prior.insert(make_pair(rSym,pow(.5,double(tree->size()))));
            //prior.insert(make_pair(rSym,.00000001));
            delete tree;
        }
        for(I2Dmap::iterator iter = tagP.begin();iter != tagP.end();++iter) {
            unsigned int rSym = iter->first;
            ParseTree* tree = getTree(rSym);
            prior.insert(make_pair(rSym,pow(.5,double(tree->size()))));
            //prior.insert(make_pair(rSym,.00000001));
            delete tree;
        }
        clearG();
    }
    
    void clearG() {
        gKids.clear();
        tKids.clear();
    }
    
    void getG() {
        clearG();
        gKids.resize(goodmanIndex);
        tKids.resize(goodmanIndex);
        for(Umap::iterator iter = umap.begin();iter != umap.end();++iter) {
            unsigned int to = iter->first;
            vector<unsigned int> kv;
            kv.push_back(to);
            set<unsigned int> froms = iter->second;
            for(set<unsigned int>::iterator jter = froms.begin();jter != froms.end();++jter) {
                unsigned int from = *jter;
                gKids[from] = kv;
            }
        }
        for(Bmap::iterator iter = bmap.begin();iter != bmap.end();++iter) {
            pair<unsigned int,unsigned int> to = iter->first;
            vector<unsigned int> kv;
            kv.push_back(to.first);
            kv.push_back(to.second);
            set<unsigned int> froms = iter->second;
            for(set<unsigned int>::iterator jter = froms.begin();jter != froms.end();++jter) {
                unsigned int from = *jter;
                gKids[from] = kv;
            }
        }
        for(S2Setmap::iterator iter = preterms.begin();iter != preterms.end();++iter) {
            string t = iter->first;
            set<unsigned int> froms = iter->second;
            for(set<unsigned int>::iterator jter = froms.begin();jter != froms.end();++jter) {
                unsigned int from = *jter;
                tKids[from] = t;
            }
        }
    }


    
    void dump() {

        printf("SYMBOLS\n");
        for(unsigned int i=0;i<goodmanIndex;++i) {
            if(i < nSym) {
                printf("%u\t: %s\n",i,syms[i].c_str());
            } else if(i < nSym*2) {
                printf("%u\t: %s*\n",i,syms[i-nSym].c_str());
            } else if(i == nSym*2) {
                printf("%u\t: %s\n",i,syms[nSym*2].c_str());
            } else {
                unsigned int bSym = baseSym[i];
                printf("%u\t: %s\n",i,syms[bSym].c_str());
            }
        }

        printf("PRETERMS\n");
        for(S2Setmap::iterator iter = preterms.begin();iter != preterms.end();++iter) {
            set<unsigned int> s = iter->second;
            for(set<unsigned int>::iterator jter = s.begin();jter != s.end();++jter) {
                unsigned int i = *jter;
                string sym = syms[baseSym[i]];
                printf("%u -> %s\t\t\t%s -> %s\n",i,iter->first.c_str(),sym.c_str(),iter->first.c_str());
            }
        }
        
        printf("UNARIES\n");
        for(Umap::iterator iter = umap.begin();iter != umap.end();++iter) {
            set<unsigned int> s = iter->second;
            for(set<unsigned int>::iterator jter = s.begin();jter != s.end();++jter) {
                unsigned int i = *jter;
                string sym = syms[baseSym[i]];
                unsigned int j = iter->first;
                string jsym = syms[baseSym[j]];
                printf("%u -> %u\t\t\t%s -> %s\n",i,j,sym.c_str(),jsym.c_str());
            }
        }

        printf("BINARIES\n");
        for(Bmap::iterator iter = bmap.begin();iter != bmap.end();++iter) {
            set<unsigned int> s = iter->second;
            for(set<unsigned int>::iterator jter = s.begin();jter != s.end();++jter) {

                unsigned int i = *jter;
                string sym = syms[baseSym[i]];
                unsigned int j = iter->first.first;
                string jsym = syms[baseSym[j]];
                unsigned int k = iter->first.second;
                string ksym = syms[baseSym[k]];

                printf("%u -> %u %u\t\t\t%s -> %s %s\n",i,j,k,sym.c_str(),jsym.c_str(),ksym.c_str());
            }
        }
    }

    //sort rules by symbol, then by descending probability
    struct Rulesort {
        
        Rulesort(BucketGrammar* gr_) : gr(gr_) {}
        
        bool operator()(pair<unsigned int,double> const &a, pair<unsigned int,double> const &b) {
            unsigned int bA = gr->baseSym[a.first];
            unsigned int bB = gr->baseSym[b.first];
            if(bA != bB)
                return bA < bB;
            else
                return a.second > b.second;
        }
        
        BucketGrammar* gr;
    };

    //sort adjunction keys by descending probability
    struct Psort {
        bool operator()(pair<string,double> const &a, pair<string,double> const &b) {
            return a.second > b.second;
        }
    };

    void trimNsave(std::ofstream& ofs, double tsgCut, double tagCut, double adjCut) {

        ofs << int(nSym) << "\n";
        for(size_t i=0;i<nSym;++i) {
            ofs << syms[i] << "\n";
        }
        
        getG();        
        
        vector<pair<unsigned int,double> > srt;
        for(I2Dmap::iterator iter = tsgP.begin();iter != tsgP.end();++iter) {
            ParseTree* tree = getTree(iter->first);
            if(iter->second > tsgCut || tree->isPCFG()) {
                srt.push_back(*iter);
            } else {
                //printf("ELIM %s\n",toString(tree).c_str());
            }
            delete tree;
        }
        sort(srt.begin(),srt.end(),Rulesort(this));
        ofs << int(srt.size()) << "\n";
        printf("%d -> %d TSG rules\n",int(nTSG),int(srt.size()));
        
        for(vector<pair<unsigned int,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            unsigned int rSym = iter->first;
            double prob = iter->second;
            ParseTree* tree = getTree(rSym);
            ofs << toString(tree) << "\n" << prob << "\n";
            delete tree;
        }

        srt.clear();

        set<unsigned int> allowB;
        for(unsigned int i =nSym*2+1;i<goodmanIndex;++i) {
            unsigned int bS = baseSym[i];
            if(ptsyms.count(bS) == 0 && bS != nSym*2 && tagP.count(i) == 0) { 
                string s = getAKey(i);
                if(aKeys.count(s) > 0 && adjoinP[aKeys[s]] > adjCut) {
                    allowB.insert(baseSym[i]);
                }
            }
        }
                
        for(I2Dmap::iterator iter = tagP.begin();iter != tagP.end();++iter) {
            ParseTree* tree = getTree(iter->first);
            if(iter->second > tagCut && allowB.count(baseSym[iter->first]) > 0) 
                srt.push_back(*iter);
            delete tree;
        }
        sort(srt.begin(),srt.end(),Rulesort(this));
        ofs << int(srt.size()) << "\n";
        printf("%d -> %d TAG rules\n",int(nTAG),int(srt.size()));
        
        for(vector<pair<unsigned int,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            unsigned int rSym = iter->first;
            double prob = iter->second;
            ParseTree* tree = getTree(rSym);
            ofs << toString(tree) << "\n" << prob << "\n";
            delete tree;
        }

        
        vector<pair<string,double> > aaa;
        for(S2Imap::iterator i = aKeys.begin(); i != aKeys.end();++i) {
            double p = adjoinP[i->second];
            if(p > adjCut)
                aaa.push_back(make_pair(i->first,p));
        }
        printf("%d -> %d Adjunction groups\n",int(aKeys.size()),int(aaa.size()));
        ofs << int(aaa.size()) << "\n";
        sort(aaa.begin(),aaa.end(),Psort());
        for(size_t i=0;i<aaa.size();++i) {
            ofs << aaa[i].first << "\n" << aaa[i].second << "\n";
        }
        
    }
    
    void write(std::ofstream& ofs) {
        
        ofs << int(nSym) << "\n";
        for(size_t i=0;i<nSym;++i) {
            ofs << syms[i] << "\n";
        }

        getG();        
        
        vector<pair<unsigned int,double> > srt;
        for(I2Dmap::iterator iter = tsgP.begin();iter != tsgP.end();++iter) {
            if(iter->second > 0)
                srt.push_back(*iter);
        }
        ofs << int(srt.size()) << "\n";
        sort(srt.begin(),srt.end(),Rulesort(this));        
        for(vector<pair<unsigned int,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            unsigned int rSym = iter->first;
            double prob = iter->second;
            ParseTree* tree = getTree(rSym);
            ofs << toString(tree) << "\n" << prob << "\n";
            delete tree;
        }

        srt.clear();
        for(I2Dmap::iterator iter = tagP.begin();iter != tagP.end();++iter) {
            if(iter->second > 0)
                srt.push_back(*iter);
        }
        ofs << int(srt.size()) << "\n";
        sort(srt.begin(),srt.end(),Rulesort(this));        
        for(vector<pair<unsigned int,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            unsigned int rSym = iter->first;
            double prob = iter->second;
            ParseTree* tree = getTree(rSym);
            ofs << toString(tree) << "\n" << prob << "\n";
            delete tree;
        }

        ofs << int(aKeys.size()) << "\n";
        vector<pair<string,double> > aaa;
        for(S2Imap::iterator i = aKeys.begin(); i != aKeys.end();++i) {            
            aaa.push_back(make_pair(i->first,adjoinP[i->second]));
        }        
        sort(aaa.begin(),aaa.end(),Psort());
        for(size_t i=0;i<aaa.size();++i) {
            ofs << aaa[i].first << "\n" << aaa[i].second << "\n";
        }
        
        clearG();
    }

    string toString(ParseTree* tree) {
        vector<string> terms = tree->terms;
        reverse(terms.begin(),terms.end());
        return toString(tree->root,terms);
    }
    
    string toString(TreeNode* node, vector<string>& terms) {
        if(baseSym[node->sym] == nSym*2) { //GLUE!!
            string ret = "";
            for(size_t i=0;i<node->k.size();++i) {
                TreeNode* k = node->k[i];
                ret += toString(k,terms);
                if(i < node->k.size()-1)
                    ret += " ";
            }
            return ret;
        } else {
            string ret = "(" + syms[baseSym[node->sym]] + " ";
            if(node->k.size() == 0) {
                ret += terms.back();
                terms.pop_back();
            } else {
                for(size_t i=0;i<node->k.size();++i) {
                    TreeNode* k = node->k[i];
                    ret += toString(k,terms);
                    if(i < node->k.size()-1)
                        ret += " ";
                }
            }
            ret += ")";
            return ret;
        }
    }

    TreeNode* recTree(unsigned int s, vector<string>& terms) {
        if(baseSym[s] == s) {
            terms.push_back("<>");
            return new TreeNode(s,vector<TreeNode*>());
        }
        vector<unsigned int> gK = gKids[s];
        if(gK.size() == 0) {
            terms.push_back(tKids[s]);
            return new TreeNode(s,vector<TreeNode*>());
        } else {
            vector<TreeNode*> k;
            for(size_t i=0;i<gK.size();++i) {
                k.push_back(recTree(gK[i],terms));
            }
            return new TreeNode(s,k);
        }
    }
    
    ParseTree* getTree(unsigned int s) {
        vector<string> terms;
        TreeNode* root = recTree(s,terms);
        return new ParseTree(root,terms);
    }

    double adjoinProb(ChartItem* it,unsigned char nBuckets) {
        if(it->buckets.size() < nBuckets) {
            return adjoinProb(it->sym);
        } else
            return 0.0;

    }

    double adjoinProb(unsigned int s) {
        //allow adjunction if 1) sym is not a base sym, 2) under limit, and 3) adjP in grammar != 0
        I2Imap::iterator fter = aClassMap.find(s);
        if(fter != aClassMap.end())
            return adjoinP[fter->second];
        else
            return 0.0;
    }
    
private:
    
};


//Accepts a set of terminals and excludes any rule with a terminal not in this set
class LimitedBucketGrammar : public BucketGrammar {
public:

    //get the words from a file
    LimitedBucketGrammar(std::ifstream& ifs, const char* allWords) : BucketGrammar() {
        std::ifstream wz(allWords);
        if(!wz.is_open()) {
            printf("Invalid file at %s\n",allWords);
            exit(-2);
        }

        dict.set_empty_key("");

        string w;
        while(wz.good()) {
            wz >> w;
            dict.insert(w);
        }
        loadGrammar(ifs);
    }

    //get the words from a space delimited string
    LimitedBucketGrammar(std::ifstream& ifs, string s) : BucketGrammar() {
        
        stringstream wz(s);

        dict.set_empty_key("");

        string w;
        while(wz.good()) {
            wz >> w;
            dict.insert(w);
        }
        
        loadGrammar(ifs);
    }
    
    bool allow(string& s) {
        bool a = true;
        vector<string> terms;
        stringstream ss(s);
        string w;
        while(ss.good()) {
            ss >> w;
            if(w[0] != '(') {

                string tok = w.substr(0,w.find_first_of(')'));
                if(tok != "<>") {

                    terms.push_back(tok);
                }
                
            }
        }
        
        for(vector<string>::iterator iter=terms.begin();a && iter!=terms.end();++iter) {
            if(dict.find(*iter) == dict.end())
                a = false;
        }
        
        return a;
    }

    google::dense_hash_set<string> dict;
    
};

#endif
