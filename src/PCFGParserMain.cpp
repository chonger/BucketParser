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


            ParseTree goldTree(goldS,parser.sym2base,parser.nSym);
            std::vector<EvalItem> gold = goldTree.getEItems(goldTree.root,parser.syms).first;
            
            PChart* chart = parser.getChart(terms);
            
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
            



            
            delete chart;

        }
    }

    double prec = tp/(tp + fp);
    double rec = tp/(tp + fn);
    printf("P : %f R : %f\n",prec,rec);
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


