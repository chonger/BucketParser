#include "BucketParser.hpp"
#include "BucketGrammar.hpp"
#include "BucketEM.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
#include "ec/Parser.h"

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
    
    if(!(argc == 5 || argc == 6)) {
        printf("parse [grammar] [nBuckets] [inputFile] [outputFile]\n");
        exit(-2);
    }

    printf("\n\n\t\tWELCOME TO THE BUCKET PARSER\n\n");
    
    const char* grammarFilename = argv[1];
    unsigned int nBuckets = atoi(argv[2]);
    const char* toparseFilename = argv[3];
    const char* outputFilename = argv[4];


    bool ctf = false;
    PCFGParser* pparser = NULL;
    if(argc == 6) {
        const char* pcfgFilename = argv[5];
        ctf = true;
    
        std::ifstream pcfgifs(pcfgFilename);
        if(!pcfgifs.is_open()) {
            printf("Invalid file at %s\n",pcfgFilename);
            exit(-2);
        }

        pparser = new PCFGParser(pcfgifs);
    }
    
    std::ifstream pifs(toparseFilename);
    if(!pifs.is_open()) {
        printf("Invalid file at %s\n",toparseFilename);
        exit(-2);
    }
    
    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }
    
    //Parser::init("/home/chonger/code/BucketParser/src/ec/DATA/");
    
    int c = 0;
    while(pifs.good() && cont) {
        //cont = false;
        q = false;
        c++;
        std::string toParse;
        getline(pifs,toParse);
        if(toParse.size() > 0) {
            std::ifstream ifs(grammarFilename);
            if(!ifs.is_open()) {
                printf("Invalid file at %s\n",grammarFilename);
                exit(-2);
            }

            //toParse = "The futures halt was even assailed by Big Board floor traders .";
            //toParse = "I ate lunch .";

            printf("[%d]\n",c); 
            
            LimitedBucketGrammar g(ifs,toParse);
            
            BucketParser* p = NULL;
            if(ctf)
                p = new CTFBucketParser(g,nBuckets,pparser);
            else
                p = new BucketParser(g,nBuckets);
            //printf("%d PARSING : %s\n",c,toParse.c_str());
            string parsed = p->parse(toParse);
            ofs << parsed << "\n";
            //printf("PARSED  : %s\n",parsed.c_str());

            
            fflush(stdout);
            
            delete p;

        }
    }

    if(ctf)
        delete pparser;
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


