#include "PCFGParser.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
using namespace std;

bool cont = true;
bool q = false;

void abortFunc(int sig) {
    if(q) 
        exit(-2);
    else {
        cont = false;
        q = true;
    }
}

struct EvalItem {

    std::string sym;
    unsigned int i,j;
                    
    EvalItem(std::string sym_, unsigned int i_, unsigned int j_) : sym(sym_), i(i_), j(j_) {
        
    }

    bool operator==(EvalItem& o) {
        return sym == o.sym && i == o.i && j == o.j;
    }
        
    
};

std::pair<std::vector<EvalItem>,unsigned int> getEItems(TreeNode* node, string* syms) {

    std::vector<EvalItem> ret;

    if(node->k.size() == 0) {
        EvalItem e(syms[node->sym],0,0);
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

        EvalItem e(syms[node->sym],tot-1,0);
        ret.push_back(e);
        return make_pair(ret,tot);
    }
    
}

int main(int argc, const char* argv[]) {

    signal(SIGABRT,&abortFunc);
    signal(SIGTERM,&abortFunc);
    signal(SIGINT,&abortFunc);
    
    printf("\n\n\t\tWELCOME TO THE BUCKET PCFG PARSER\n\n");
    
    const char* grammarFilename = argv[1];
    const char* toparseFilename = argv[2];
    const char* goldFilename = argv[3];

    std::ifstream pifs(grammarFilename);
    if(!pifs.is_open()) {
        printf("Invalid file at %s\n",grammarFilename);
        exit(-2);
    }

    PCFGParser parser(pifs);

    std::ifstream ifs(toparseFilename);
    if(!ifs.is_open()) {
        printf("Invalid file at %s\n",toparseFilename);
        exit(-2);
    }

    std::ifstream gifs(goldFilename);
    if(!gifs.is_open()) {
        printf("Invalid file at %s\n",goldFilename);
        exit(-2);
    }

    double tp = 0.0;
    double fp = 0.0;
    double fn = 0.0;
    
    
    std::string toParse = "";
    std::string goldS = "";
    while(ifs.good() && cont) {
        q = false;
        getline(ifs,toParse);
        getline(gifs,goldS);
        if(toParse.size() > 0) {
            
            vector<string> terms;
            stringstream ss;
            ss << toParse;
            while(ss.good()) {
                string t;
                ss >> t;
                terms.push_back(t);
            }

            //            printf("%s\n",goldS.c_str());
            ParseTree goldTree(goldS,parser.sym2base,parser.nSym);
            std::vector<EvalItem> gold = getEItems(goldTree.root,parser.syms).first;
            
            PChart* chart = parser.getChart(terms);
            
            std::vector<EvalItem> predict;
            
            for(int i=0;i<chart->size;++i) {
                int rowS = chart->size-i;
                for(int j=0;j<rowS;++j) {
                    PCell& pc = chart->cells[i][j];
                    if(!pc.empty()) {
                        for(PCell::iterator iter = pc.begin();iter != pc.end();++iter) {
                            PData* d = iter->second;
                            if(d->prob >= 0.001) {
                                EvalItem e(parser.syms[d->sym],i,j);
                                predict.push_back(e);
                            }
                        }
                    }
                }
            }
            /**
            printf("GOLD\n");
            for(size_t i=0;i<gold.size();++i) {
                EvalItem& e = gold[i];
                printf("%d %d %s\n",e.i,e.j,e.sym.c_str());
            }
            
            printf("PREDICT\n");
            for(size_t i=0;i<predict.size();++i) {
                EvalItem& e = predict[i];
                printf("%d %d %s\n",e.i,e.j,e.sym.c_str());
            }
            */
            double mytp = 0.0;
            for(size_t i=0;i<gold.size();++i) {
                EvalItem& e1 = gold[i];
                for(size_t j=0;j<predict.size();++j) {
                    EvalItem& e2 = predict[j];
                    if(e1 == e2)
                        mytp += 1.0;
                }
            }

            double myfp = predict.size() - mytp;
            double myfn = gold.size() - mytp;

            tp += mytp;
            fp += myfp;
            fn += myfn;
            
            double prec = tp/(tp + fp);
            double rec = tp/(tp + fn);


            
            delete chart;

        }
    }

    printf("P : %f R : %f\n",prec,rec);
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


