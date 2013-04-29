#include "BucketVB.hpp"
#include "VBGrammar.hpp"
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

    srand(12345);
    
    printf("\n\n\t\tVARIATIONAL INFERENCE FOR TSGs\n\n");
    
    const char* grammarFilename = argv[1];
    const char* trainFilename = argv[2];
    unsigned int nIter = atoi(argv[3]);
    const char* outputFilename = argv[4];
    
    std::ifstream ifs(grammarFilename);
    if(!ifs.is_open()) {
        printf("Invalid file at %s\n",grammarFilename);
        exit(-2);
    }

    std::ifstream tifs(trainFilename);
    if(!tifs.is_open()) {
        printf("Invalid file at %s\n",trainFilename);
        exit(-2);
    }

    printf("%u iterations\n",nIter);
    
    VBGrammar g(ifs,tifs);

    printf("Grammar made\n");
    
    BucketVB vb(g);
    
    
    vb.emIter(true);
    vb.reformGrammar();
    vb.emIter(true);
    
    /**
    for(unsigned int i=0;i<100 && cont;++i) {
            printf("%u - %e\n",i+1,vb.emIter(true));
    }
    */
    /**
    for(unsigned int i=0;i<2 && cont;++i) {
        for(unsigned int j=0;j<2;++j) {
            q = false;
            printf("%u - %e\n",i+1,vb.emIter(true));
        }
        vb.reformGrammar();
    }
     
    for(unsigned int j=0;j<2;++j) {
        q = false;
        printf("%e\n",vb.emIter(true));
    }
    */
    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }
    g.writeTruncated(ofs);
    ofs.close();

    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


