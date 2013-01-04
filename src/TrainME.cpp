#include "PCFGParser.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
#include <maxent/maxentmodel.hpp>
using namespace std;
using namespace maxent;

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
    
    MaxentModel m;
    verbose = 1;
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

    string goldS = "";

    m.begin_add_event();
    int i=0;
    while(ifs2.good() && cont) {
        q = false;
        ++i;
        printf("[%d]\n",i);
        getline(ifs2,goldS);
        ParseTree goldTree(goldS,parser.sym2base,parser.nSym);
        vector<EvalItem> gold = goldTree.getEItems(goldTree.root,parser.syms).first;
        PChart* chart = parser.getChart(goldTree.terms);
        vector<EvalItem> predict = chart->getEItems(parser.syms);

        for(size_t i=0;i<gold.size();++i) {
            EvalItem& e = gold[i];
            vector<pair<string,float> > context;
            for(size_t j=0;j<predict.size();++j) {
                EvalItem& p = predict[j];
                if(p.i == e.i && p.j == e.j) {
                    context.push_back(make_pair(p.sym,float(p.prob)));
                }
            }
            m.add_event(context,e.sym);
        }
    }
        
    m.end_add_event();

    m.train();

    m.save(outputFilename);
    
}
