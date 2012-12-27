#include "BucketEM.hpp"
#include "BucketGrammar.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>

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

    BucketGrammar g(ifs);
    
    std::ofstream ofs(outputFilename);
    if(!ofs.is_open()) {
        printf("Invalid file at %s\n",outputFilename);
        exit(-2);
    }
    
    g.trimNsave(ofs,.000000001,.00000001,.001);
    
    printf("\n\n\t\tTHANK YOU FOR PLAYING\n\n\t\t\tTHE END\n\n");
    
    return 0;
    
}


