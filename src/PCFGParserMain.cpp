#include "PCFGParser.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
#include <math.h>

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

    
    size_t nTrials = 10;
    double tp[nTrials];
    double fp[nTrials];
    double fn[nTrials];    
    double cuts[nTrials];
    
    for(size_t i=0;i<nTrials;++i) {
        tp[i] = 0.0;
        fp[i] = 0.0;
        fn[i] = 0.0;
        double exp = i;
        exp += 1.0;
        exp *= -1;
        cuts[i] = pow(10.0,exp);
    }

    
    std::string toParse = "";
    std::string goldS = "";
    int iii = 0;
    while(ifs.good() && cont) {
        q = false;
        getline(ifs,toParse);
        getline(gifs,goldS);

        iii += 1;
        printf("[%d]\n",iii);
        if(toParse.size() > 0) {
            
            vector<string> terms;
            stringstream ss;
            ss << toParse;
            while(ss.good()) {
                string t;
                ss >> t;
                terms.push_back(t);
            }
            
            ParseTree goldTree(goldS,parser.sym2base,parser.nSym);
            
            std::vector<EvalItem> gold = goldTree.getEItems(goldTree.root,parser.syms).first;
            
            PChart* chart = parser.getChart(terms);

            vector<EvalItem> predict = chart->getEItems(parser.syms);
            
            EvalMap em;
            em.set_empty_key(make_pair(-1,-1));
            for(size_t j=0;j<predict.size();++j) {
                EvalItem& p = predict[j];
                em[make_pair(p.i,p.j)].push_back(&p);
            }

            EvalMap emG;
            emG.set_empty_key(make_pair(-1,-1));
            for(size_t j=0;j<gold.size();++j) {
                EvalItem& p = gold[j];
                emG[make_pair(p.i,p.j)].push_back(&p);
            }

            /**            
            std::vector<EvalItem> predict;
            for(int i=0;i<chart->size;++i) {
                int rowS = chart->size-i;
                for(int j=0;j<rowS;++j) {
                    PCell& pc = chart->cells[i][j];
                    if(!pc.empty()) {
                        for(PCell::iterator iter = pc.begin();iter != pc.end();++iter) {
                            PData* d = iter->second;
                            if(d->prob >= .0001) {
                                EvalItem e(parser.syms[d->sym],i,j,d->prob);
                                predict.push_back(e);
                            }
                        }
                    }
                }
            }
            */
            /**            
            printf("%s\n",goldS.c_str());
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
            double mytp[nTrials];
            double myfp[nTrials];
            double myfn[nTrials];

            for(size_t i=0;i<nTrials;++i) {
                mytp[i] = 0.0;
                myfp[i] = 0.0;
                myfn[i] = 0.0;
            }

            for(EvalMap::iterator iter = emG.begin();iter != emG.end();++iter) {
                pair<unsigned int,unsigned int> k = iter->first;
                vector<EvalItem*>& gs = iter->second;
                vector<EvalItem*>& ps = em[k];

                vector<EvalItem*> preds[nTrials];
                for(size_t i=0;i<nTrials;++i) {
                    for(size_t j=0;j<ps.size();++j) {
                        if(ps[j]->prob > cuts[i]) {
                            preds[i].push_back(ps[j]);
                        }
                    }
                }
                for(size_t j=0;j<nTrials;++j) {
                    vector<EvalItem*>& ppp = preds[j];
                    double found = 0.0;

                    for(size_t i=0;i<gs.size();++i) {
                        EvalItem* e = gs[i];
                        for(size_t k=0;k<ppp.size();++k) {
                            if(ppp[k]->sym == e->sym) {
                                found += 1.0;
                                break;
                            }
                        }
                    }

                    mytp[j] += found;
                    myfp[j] += ppp.size() - found;
                }
            }

            /**
            for(size_t i=0;i<gold.size();++i) {
                EvalItem& e = gold[i];
                vector<EvalItem*>& ps = em[make_pair(e.i,e.j)];
                for(size_t j=0;j<ps.size();++j) {
                    EvalItem* p = ps[j];
                    for(size_t i=0;i<nTrials;++i) {
                        if(p->prob > cuts[i]) {
                            if(p->sym == e.sym)
                                mytp[i] += 1.0;
                            else
                                myfp[i] += 1.0;
                        }
                    }
                }
            }
            */
            for(size_t i=0;i<nTrials;++i) {
                double myfn = gold.size() - mytp[i];
                tp[i] += mytp[i];
                fp[i] += myfp[i];
                fn[i] += myfn;
            }
            
            delete chart;
            
            for(size_t i=0;i<nTrials;++i) {
                double prec = tp[i]/(tp[i] + fp[i]);
                double rec = tp[i]/(tp[i] + fn[i]);
                printf("%E %f %f\n",cuts[i],prec,rec);
            }
            printf("\n");
        }
    }

    for(size_t i=0;i<nTrials;++i) {
        double prec = tp[i]/(tp[i] + fp[i]);
        double rec = tp[i]/(tp[i] + fn[i]);
        printf("%E %f %f\n",cuts[i],prec,rec);
    }
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


