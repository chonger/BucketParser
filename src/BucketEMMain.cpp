#include "BucketEM.hpp"
#include "BucketGrammar.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>

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
    
    if(argc != 6) {
        printf("parse [grammar] [nBuckets] [nIter] [inputFile] [outputFile]\n");
        exit(-2);
    }

    printf("\n\n\t\tWELCOME TO THE BUCKET PARSER\n\n");
    
    const char* grammarFilename = argv[1];
    unsigned int nBuckets = atoi(argv[2]);
    unsigned int nIter = atoi(argv[3]);
    const char* trainFilename = argv[4];
    const char* outputFilename = argv[5];
    std::ifstream ifs(grammarFilename);
    if(!ifs.is_open()) {
        printf("Invalid file at %s\n",grammarFilename);
        exit(-2);
    }

    printf("%u buckets, %u iterations of EM\n",nBuckets,nIter);
    
    BucketGrammar g(ifs);
    BucketEM em(g,trainFilename,nBuckets);

    for(unsigned int i=0;i<nIter && cont;++i) {
        q = false;
        printf("%u - %e\n",i+1,em.emIter(true));
    }

    em.emIter(false);
    
    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }

    g.write(ofs);
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


