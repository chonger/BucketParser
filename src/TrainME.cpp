#include "PCFGParser.hpp"
#include <iostream>
#include <stdio.h>
#include <string>
#include <climits>
#include <signal.h>
#include <maxent/maxentmodel.hpp>
using namespace std;
using namespace maxent;

int main(int argc, const char* argv[]) {

    MaxentModel m;

    m.begin_add_event();

    vector<string> context1;
    context1.push_back("A");
    context1.push_back("C");

    vector<string> context2;
    context2.push_back("B");
    context2.push_back("C");

    m.add_event(context1,"YES");
    m.add_event(context2,"NO");
    
    m.end_add_event();

    m.train();

    printf("%f\n",m.eval(context1,"YES"));
    printf("%f\n",m.eval(context2,"YES"));
    
    printf("DONE\n");
    
}
