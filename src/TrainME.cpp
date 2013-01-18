#include "PCFGParser.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
//#include <maxent/maxentmodel.hpp>
using namespace std;
//using namespace maxent;

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
    

    const char* pcfgFilename = argv[1];
    const char* trainFilename = argv[2];
    const char* outputFilename = argv[3];


    
    std::ifstream ifs1(pcfgFilename);
    if(!ifs1.is_open()) {
        printf("Invalid file at %s\n",pcfgFilename);
        exit(-2);
    }
    PCFGParser parser(ifs1);

    std::ifstream ifs2(trainFilename);
    if(!ifs2.is_open()) {
        printf("Invalid file at %s\n",trainFilename);
        exit(-2);
    }

    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }
    
    size_t nS = parser.nSym;
    //MaxentModel m[nS];
    //  verbose = 1;
    
    string goldS = "";
    /**
    for(size_t i=0;i<nS;++i) {
        m[i].begin_add_event();
    }
    */
    int i=0;
    while(ifs2.good() && cont) {
        q = false;
        ++i;
        printf("[%d]\n",i);
        getline(ifs2,goldS);
        
        //        printf("%s\n",goldS.c_str());
        ParseTree goldTree(goldS,parser.sym2base,parser.nSym);
        vector<EvalItem> gold = goldTree.getEItems(goldTree.root,parser.syms).first;
        PChart* chart = parser.getChart(goldTree.terms);
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
        
        for(EvalMap::iterator iter = em.begin();iter != em.end();++iter) {
            
            pair<unsigned int,unsigned int> ind = iter->first;
            vector<EvalItem*>& ps = iter->second;

            //printf("\n%d %d GOLD : ",ind.first,ind.second);
            set<string> golds;
            vector<EvalItem*> gs = emG[ind];
            for(size_t i=0;i<gs.size();++i) {
                golds.insert(gs[i]->sym);
                //printf("%s ",gs[i]->sym.c_str());
                ofs << gs[i]->sym << " ";
            }
            //printf("\n");
            ofs << ": ";
            vector<pair<string,float> > context;
            
            for(size_t j=0;j<ps.size();++j) {
                EvalItem* p = ps[j];
                context.push_back(make_pair(p->sym,float(p->prob)));
                ofs << p->sym << " " << p->prob << " ";
                /**
                if(p->prob > .0001) {
                    printf("%s %f\n",p->sym.c_str(),p->prob);
                }
                */
            }

            ofs << "\n";
            

            /**
            for(unsigned int i=0;i<nS;++i) {
                if(golds.count(parser.syms[i]) > 0) 
                    m[i].add_event(context,"YES");
                else
                    m[i].add_event(context,"NO");
            }
            */

        }

        
        delete chart;
    }
    /**    
    for(size_t i=0;i<nS;++i) {        
        m[i].end_add_event();
        m[i].train();
        string o = outputFilename + parser.syms[i];
        m[i].save(o);
    }
    */
    
}
