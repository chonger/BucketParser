#ifndef VB_GRAMMAR
#define VB_GRAMMAR

#import "BucketGrammar.hpp"

/**
 *   Strategy - 
 *
 *     1) input a format that is the base PCFG (standard format)
 *        followed by a list of current TSG rules with their SB weights,
 *        this also defines the mass for "the rest" for any symbol
 *        
 *
 */ 

typedef google::dense_hash_map<vector<unsigned int>,double,hash<vector <unsigned int> >, equal_to<vector <unsigned int> > > V2Dmap;
typedef google::dense_hash_map<pair<string, unsigned int>,double,hash<pair<string , unsigned int> >,equal_to<pair<string, unsigned int> > > SIDmap;
typedef google::dense_hash_map<vector<unsigned int>,unsigned int,hash<vector <unsigned int> >, equal_to<vector <unsigned int> > > V2Pmap;
typedef google::dense_hash_map<pair<string, unsigned int>,unsigned int,hash<pair<string , unsigned int> >,equal_to<pair<string, unsigned int> > > SIPmap;

typedef google::dense_hash_map<ParseTree*,unsigned int,ParseTreeHash,ParseTreeEq > T2Imap;
typedef google::dense_hash_map<ParseTree*,double,ParseTreeHash,ParseTreeEq > T2Dmap;

class VBGrammar {
public:

    VBGrammar(std::ifstream& ifs,std::ifstream& trainifs) {
        
        //HASMAP SETUP
        //these hashmaps are data members, and require this setup just once
        tsgP.set_empty_key(-1);
        prior.set_empty_key(-1);

        vector<unsigned int> vectorkey1;
        vector<unsigned int> vectorkey2;
        vectorkey2.push_back(1);
        
        rulemap.set_empty_key(vectorkey1);
        rulemap.set_deleted_key(vectorkey2);

        eTree = new ParseTree("EMPTY");
        dTree = new ParseTree("DELETED");
        ruleTreez.set_empty_key(eTree);
        ruleTreez.set_deleted_key(dTree);
        
        preterms.set_empty_key(make_pair("",0));
        preterms.set_deleted_key(make_pair("",1));

        sym2base.set_empty_key("EMPTYKEY");
        sym2base.set_deleted_key("DELETEDKEY");

        pcfgPT.set_empty_key(make_pair("",0));
        pcfgNT.set_empty_key(vector<unsigned int>());
        
        //Start reading the grammar file here
        string line;
        getline(ifs,line);
        nSym = strtol(line.c_str(),NULL,0);
        printf("%d syms\n",nSym);
        syms = new string[nSym];
        cutProb = new double[nSym];
        alpha = new double[nSym];
        
        for(int i=0;i<nSym;++i) {
            getline(ifs,line);
            string sym = line;
            if(sym == "ROOT")
                rootSym = i;
            syms[i] = sym;
            sym2base[sym] = i;
            cutProb[i] = .5;
        }
        for(unsigned int i=0;i<nSym;++i) { //map the base syms to themselves
            baseSym.push_back(i);
        }
        for(unsigned int i=nSym;i<nSym*2;++i) { //map the primed syms to the base
            baseSym.push_back(i-nSym);
        }

        for(size_t i=0;i<nSym;++i) {
            printf("%s\n",syms[i].c_str());
        }

        printf("Reading training trees\n");
        while(trainifs.good()){
            string treeS;
            getline(trainifs,treeS);
            if(treeS.size() > 0) {
                //printf("TREE: %s\n",treeS.c_str());
                trainTreez.push_back(growTree(treeS));
            }
        }
        trainifs.close();
        
        //build prime grammar...
        loadPrime();

        //read TSG rules        
        getline(ifs,line);
        nTSG = strtol(line.c_str(),NULL,0);
        printf("Initializing with %d TSG rules\n",nTSG);

        vector<ParseTree*> initTreez;
        vector<double> initProbs;
        
        for(int i=0;i<nTSG;++i) {
            getline(ifs,line); 
            string s = line;
            ParseTree* pt = growTree(s);
            initTreez.push_back(pt);
            getline(ifs,line);
            double prob = atof(line.c_str());
            assert(prob <= 1.0 && prob >= 0.0);
            initProbs.push_back(prob);
        }
        ifs.close();

        makeTransform(initTreez,initProbs);

        /**
        for(size_t i=0;i<initTreez.size();++i) {
            delete initTreez[i];
        }
        */
        
        printf("done constructing\n");
    }

    ~VBGrammar() {

        for(vector<ParseTree*>::iterator iter = trainTreez.begin();iter != trainTreez.end();++iter) {
            delete *iter;
        }

        for(T2Imap::iterator iter = ruleTreez.begin();iter != ruleTreez.end();++iter) {
            delete iter->first;
        }
        
        delete[] syms;
        delete[] cutProb;
        delete[] alpha;

        delete eTree;
        delete dTree;
    }
    
    void loadPrime() {

        //get the PCFG probs - fill pcfgPT and pcfgNT

        double totals[nSym];
        for(size_t i=0;i<nSym;++i)
            totals[i] = 0.0;

        
        printf("%d training trees\n",trainTreez.size());
        for(size_t i=0;i<trainTreez.size();++i) {
            ParseTree* tree = trainTreez[i];
            vector<pair<TreeNode*,int> > ns = tree->nodes();
            for(size_t i=0;i<ns.size();++i) {
                TreeNode* n = ns[i].first;
                totals[n->sym] += 1.0;
                if(n->k.size() == 0) {//terminal
                    string term = tree->terms[ns[i].second];
                    pcfgPT[make_pair(term,n->sym)] += 1.0;
                } else {
                    vector<unsigned int> key;
                    key.push_back(n->sym);
                    for(size_t j=0;j<n->k.size();++j) {
                        key.push_back(n->k[j]->sym);
                    }
                    pcfgNT[key] += 1.0;
                }
            }
        }

        //normalize
        for(SIDmap::iterator iter = pcfgPT.begin();iter != pcfgPT.end();++iter) {
            pair<const pair<string,unsigned int>,double>& e = *iter;
            unsigned int s = e.first.second;
            e.second /= totals[s];
            //            printf("%s %s - %f\n",e.first.first.c_str(),syms[e.first.second].c_str(),e.second);
        }
        
        for(V2Dmap::iterator iter = pcfgNT.begin();iter != pcfgNT.end();++iter) {
            pair<const vector<unsigned int>,double>& e = *iter;
            unsigned int s = e.first[0];
            e.second /= totals[s];
            /**
            for(size_t i=0;i<e.first.size();++i) {
                printf("%s ",syms[e.first[i]].c_str());
            }
            printf(" %f\n",e.second);
            */
        }
        
    }
    
        
    void makeTransform(vector<ParseTree*> tsgs, vector<double> probs) {
        
        //this one maps terminals to the base symbols
        S2Setmap pt2base;
        pt2base.set_empty_key("EMPTYKEY");
        pt2base.set_deleted_key("DELETEDKEY");

        goodmanIndex = nSym*2; //skip both base and prime nodes

        //clear old maps
        rulemap.clear();
        ruleTreez.clear();
        preterms.clear();
        sym2base.clear();
        baseSym.clear();
        tsgP.clear();
        prior.clear();
        
        for(unsigned int i=0;i<nSym;++i) { //map the base syms to themselves
            baseSym.push_back(i);
        }
        for(unsigned int i=nSym;i<nSym*2;++i) { //map the primed syms to the base
            baseSym.push_back(i-nSym);
        }

        for(size_t i=0;i<nSym;++i) {
            alpha[i] = 1.0;
        }
        
        for(size_t i=0;i<tsgs.size();++i) {
            ParseTree* pt = tsgs[i];
            int tInd = 0;
            vector<string>::iterator titer = pt->terms.begin();
            unsigned int index = getRule(pt->root,titer);
            tsgP[index] = probs[i];
            assert(probs[i] > 0);
            ruleTreez[pt] = index;
            //printf("index tsg - %u\n",index);
            alpha[baseSym[index]] -= probs[i];
            if(alpha[baseSym[index]] <= 0)
                printf("TOO MUCH PROB at %s - %f\n",syms[baseSym[index]].c_str(),alpha[baseSym[index]]);
            assert(alpha[baseSym[index]] >= 0);
            
            //calculate P_0 prob (unnecessary?)
            double p0 = getPrior(pt);
            prior[index] = p0;
        }
        /**
        printf("Transform keys\n");
        for(V2Pmap::iterator iter = rulemap.begin();iter != rulemap.end();++iter) {
            vector<unsigned int> key = iter->first;
            for(size_t i=0;i<key.size();++i)
                printf("%u ",key[i]);
            printf("\n");
        }
        */
    }

    double getPrior(ParseTree* tree) {
        double prior = 1.0;
        vector<pair<TreeNode*,int> > ns= tree->nodes();
        for(size_t i=0;i<ns.size();++i) {
            TreeNode* n = ns[i].first;
            if(n->k.size() == 0) {//terminal
                string term = tree->terms[ns[i].second];
                if(term == "<>")
                    prior *= cutProb[n->sym];
                else
                    prior *= pcfgPT[make_pair(term,n->sym)] * (1-cutProb[n->sym]);
            } else {
                vector<unsigned int> key;
                key.push_back(n->sym);
                for(size_t j=0;j<n->k.size();++j) {
                    key.push_back(n->k[j]->sym);
                }
                prior *= pcfgNT[key] * (1-cutProb[n->sym]);
            }
        }
    }
    
    vector<string> getTerms(string& s) {
        vector<string> terms;
        stringstream ss(s);
        string w;
        while(ss.good()) {
            ss >> w;
            if(w[0] != '(') {
                string tok = w.substr(0,w.find_first_of(')'));                
                terms.push_back(tok);
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
        //        printf("%s\n",sym.c_str());
        unsigned int symI = sym2base[sym];
        //printf("%u\n",symI);
        
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

    ParseTree* growTree(std::string s) {
        stringstream ss(s);
        TreeNode* root = readNode(ss);
        vector<string> terms = getTerms(s);
        ParseTree* pt = new ParseTree(root,terms);
        pt->setHashString(syms);
        return pt;
    };

    unsigned int getRule(TreeNode* node, vector<string>::iterator& titer) {

        unsigned int symI = node->sym;        
        vector<unsigned int> kSyms; //key to the transform
        kSyms.push_back(symI);
        
        if(node->k.size()==0) { //its a terminal or ntleaf
            string term = *titer;
            titer++;
            //printf("TERM %s\n",term.c_str());
            if(term == "<>") { //nonterminal leaf
                return symI;
            } else { //preterminal node
                pair<string,unsigned int> key = make_pair(term,symI);
                SIPmap::iterator fter = preterms.find(key);
                if(fter != preterms.end()) {
                    return fter->second;
                }
                //if we get here, its a new preterm
                unsigned int index = goodmanIndex;
                preterms[key] = index;
                baseSym.push_back(symI);
                ++goodmanIndex;
                return index;
            }
        } else {       
            for(size_t i=0;i<node->k.size();++i) {
                TreeNode* kNode = node->k[i];
                unsigned int kIndex = getRule(kNode,titer);
                kSyms.push_back(kIndex);                
            }

            unsigned int index = 0;
            V2Pmap::iterator fter = rulemap.find(kSyms);
            if(fter == rulemap.end()) { //new node
                index = goodmanIndex;
                rulemap[kSyms] = index;
                baseSym.push_back(symI);
                goodmanIndex++;
            } else { //seen it already
                index = fter->second;
            }
            return index;
        }
    }
    
    unsigned int primeSym(unsigned int n) {
        return baseSym[n] + nSym;
    }

    bool isPrime(unsigned int n) {
        return n >= nSym && n < 2*nSym;
    }

    struct Rulesort2 {
        
        Rulesort2() {}
        
        bool operator()(pair<ParseTree*,double> const &a, pair<ParseTree*,double> const &b) {
            unsigned int bA = a.first->root->sym;
            unsigned int bB = b.first->root->sym;
            if(bA != bB)
                return bA > bB;
            else
                return a.second > b.second;
        }
        
    };
    
    void write(std::ofstream& ofs) {
        
        ofs << int(nSym) << "\n";
        for(size_t i=0;i<nSym;++i) {
            ofs << syms[i] << "\n";
        }
        
        vector<pair<ParseTree*,double> > srt;
        for(T2Imap::iterator iter = ruleTreez.begin();iter != ruleTreez.end();++iter) {
            srt.push_back(make_pair(iter->first,tsgP[iter->second]));
        }
        ofs << int(srt.size()) << "\n";
        sort(srt.begin(),srt.end(),Rulesort2());        
        for(vector<pair<ParseTree*,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            ParseTree* t = iter->first;
            double prob = iter->second;
            ofs << t->hashv << "\n" << prob << "\n";
        }

    }
    
    V2Dmap pcfgNT;
    SIDmap pcfgPT;
    vector<ParseTree*> trainTreez;
    T2Imap ruleTreez;
    S2Imap sym2base;
    string* syms;
    unsigned int goodmanIndex;
    unsigned int rootSym;
    unsigned int nTSG;
    I2Dmap tsgP;
    I2Dmap prior;
    vector<unsigned int> baseSym;
    SIPmap preterms; //maps terminal strings to sets of goodman transformed preterms
    V2Pmap rulemap;
    unsigned int nSym;
    double* cutProb;
    double* alpha;


    ParseTree* eTree;
    ParseTree* dTree;
};


#endif
