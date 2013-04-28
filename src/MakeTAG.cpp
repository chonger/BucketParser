#include "BucketEM.hpp"
#include "BucketGrammar.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>

class TAGExtractor : public BucketGrammar {
public:
    
    TAGExtractor(std::ifstream& ifs) : BucketGrammar(ifs) {}

    void extractTAG(std::ofstream& ofs) {
        
        ofs << int(nSym) << "\n";
        for(size_t i=0;i<nSym;++i) {
            ofs << syms[i] << "\n";
        }

        getG();
        set<string> tagS,tsgS;        

        vector<pair<unsigned int,double> > srt;
        for(I2Dmap::iterator iter = tsgP.begin();iter != tsgP.end();++iter) {
            srt.push_back(*iter);
        }
        sort(srt.begin(),srt.end(),Rulesort(this));        
        for(vector<pair<unsigned int,double> >::iterator iter= srt.begin();iter!=srt.end();++iter) {
            unsigned int rSym = iter->first;
            //double prob = iter->second;
            ParseTree* tree = getTree(rSym);
            tsgS.insert(toString(tree));
            //ofs << toString(tree) << "\n" << prob << "\n";

            pair<vector<ParseTree*>,vector<ParseTree*> >  res = getTAGs(tree);

            //printf("Got %lu rules from %s\n",res.first.size(),toString(tree).c_str());

            vector<ParseTree*>& tsgs = res.first;
            vector<ParseTree*>& tags = res.second;
            
            for(size_t i=0;i<tags.size();++i) {
                //printf("TAG! - %s\n",toString(tags[i]).c_str());
                tagS.insert(toString(tags[i]));
                delete tags[i];
            }
            for(size_t i=0;i<tsgs.size();++i) {
                if(tsgs[i]->depth() > 0) {
                    //printf("TSG! %lu - %s\n",tsgs[i]->depth(),toString(tsgs[i]).c_str());
                    tsgS.insert(toString(tsgs[i]));
                }
                delete tsgs[i];
            }
            
            delete tree;
        }

        printf("%lu TSG -> %lu TSG and %lu TAG\n",srt.size(),tsgS.size(),tagS.size());
        
        ofs << int(tsgS.size()) << "\n";
        for(set<string>::iterator iter = tsgS.begin();iter != tsgS.end();++iter) {
            ofs << *iter << "\n.1\n";
        }        
        
        ofs << int(tagS.size()) << "\n";
        for(set<string>::iterator iter = tagS.begin();iter != tagS.end();++iter) {
            ofs << *iter << "\n.1\n";
        }

        getG();
        set<string> kset;
        for(unsigned int i =nSym*2+1;i<goodmanIndex;++i) {
            unsigned int bS = baseSym[i];
            if(ptsyms.count(bS) == 0 && bS != nSym*2 && tagP.count(i) == 0) { //get rid of preterminals,GLUE symbol,and tagroots
                string k = getAKey(i);
                kset.insert(k);
            }
        }                
        ofs << int(kset.size()) << "\n";
        for(set<string>::iterator i = kset.begin(); i != kset.end();++i) {
            ofs << *i << "\n.5\n";
        }
        
    }
    
    
    vector<TreeNode*> below(TreeNode* n, unsigned int s) {
        vector<TreeNode*> ret;
        if(baseSym[s] == baseSym[n->sym]) {
            ret.push_back(n);
        }
        for(vector<TreeNode*>::iterator iter = n->k.begin();iter != n->k.end();++iter) {
            vector<TreeNode*> kBelow = below(*iter,s);
            for(vector<TreeNode*>::iterator jter = kBelow.begin();jter != kBelow.end();++jter) {
                ret.push_back(*jter);
            }
        }
        return ret;
    }

    pair<vector<ParseTree*>, vector<ParseTree*> > getTAGs(ParseTree* tree) {
        vector<ParseTree*> retTAG,retTSG;
        vector<TreeNode*> ns = tree->root->nodes();
        for(size_t i=0;i<ns.size();++i) {
            TreeNode* n = ns[i];
            if(n->sym != nSym*2) {
                vector<TreeNode*> belows = below(n,n->sym);
                for(size_t j=1;j<belows.size();++j) {
                    vector<string> newtermsTAG,newtermsTSG;
                    vector<string> oldtermsTAG = tree->terms;
                    vector<string> oldtermsTSG = tree->terms;
                    reverse(oldtermsTAG.begin(),oldtermsTAG.end());
                    reverse(oldtermsTSG.begin(),oldtermsTSG.end());
                    TreeNode* newTAG = tagcopy(tree->root,n,belows[j],newtermsTAG,oldtermsTAG,false);
                    TreeNode* newTSG = tsgcopy(tree->root,n,belows[j],newtermsTSG,oldtermsTSG,true);

                    assert(newTAG != NULL);
                    assert(newTSG != NULL);

                    retTAG.push_back(new ParseTree(newTAG,newtermsTAG));
                    retTSG.push_back(new ParseTree(newTSG,newtermsTSG));
                }
            }
        }
        return make_pair(retTSG,retTAG);
    }

    TreeNode* tsgcopy(TreeNode* n, TreeNode* head, TreeNode* foot, vector<string>& newterms,vector<string>& oldterms,bool on) {
        if(on) {
            if(n == head) {
                return tsgcopy(n,head,foot,newterms,oldterms,false); 
            } else {
                if(n->k.size() == 0) {
                    newterms.push_back(oldterms.back());
                    oldterms.pop_back();
                }
                vector<TreeNode*> k;
                for(size_t i=0;i<n->k.size();++i) {
                    k.push_back(tsgcopy(n->k[i],head,foot,newterms,oldterms,true));
                }
                return new TreeNode(n->sym,k);
            }
        } else {
            if(n == foot) { //turn copying on
                return tsgcopy(n,head,foot,newterms,oldterms,true); 
            } else {
                if(n->k.size() == 0) { //burn the terminal
                    oldterms.pop_back();
                }

                TreeNode* rt = NULL;
                
                for(size_t i=0;i<n->k.size();++i) {
                    TreeNode* tn = tsgcopy(n->k[i],head,foot,newterms,oldterms,false);
                    if(tn != NULL)
                        rt = tn;
                }
                return rt;
            }
        }
    }
    
    TreeNode* tagcopy(TreeNode* n, TreeNode* head, TreeNode* foot, vector<string>& newterms,vector<string>& oldterms,bool on) {
        if(on) {
            if(n == foot) {
                newterms.push_back("<>");
                tagcopy(n,head,foot,newterms,oldterms,false); //burn through dominated terms
                return new TreeNode(baseSym[foot->sym] + nSym,vector<TreeNode*>());
            } else {
                if(n->k.size() == 0) {
                    newterms.push_back(oldterms.back());
                    oldterms.pop_back();
                }
                vector<TreeNode*> k;
                for(size_t i=0;i<n->k.size();++i) {
                    k.push_back(tagcopy(n->k[i],head,foot,newterms,oldterms,true));
                }
                return new TreeNode(n->sym,k);
            }
        } else {
            if(n == head) { //turn copying on
                return tagcopy(n,head,foot,newterms,oldterms,true); 
            } else {
                if(n->k.size() == 0) { //burn the terminal
                    oldterms.pop_back();
                }

                TreeNode* rt = NULL;
                for(size_t i=0;i<n->k.size();++i) {
                    TreeNode* tn = tagcopy(n->k[i],head,foot,newterms,oldterms,false);
                    if(tn != NULL)
                        rt = tn;
                }
                return rt;
            }
        }
    }

};

int main(int argc, const char* argv[]) {

    if(argc != 3) {
        printf("gettag [grammar] [outputFile]\n");
        exit(-2);
    }

    printf("\n\n\t\tWELCOME TO THE BUCKET PARSER\n\n");
    
    const char* grammarFilename = argv[1];
    const char* outputFilename = argv[2];
    std::ifstream ifs(grammarFilename);
    if(!ifs.is_open()) {
        printf("Invalid file at %s\n",grammarFilename);
        exit(-2);
    }

    TAGExtractor g(ifs);
    
    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }

    g.extractTAG(ofs);
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


