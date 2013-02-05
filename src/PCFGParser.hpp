#ifndef BUCKET_PCFG_PARSER
#define BUCKET_PCFG_PARSER 1

#include <cstddef>
#include <string>
#include <sstream>
#include <vector>
#include "BucketGrammar.hpp"
#include "Chart.hpp"
#include "ec/Parser.h"

typedef google::dense_hash_map<unsigned int,set<pair<unsigned int, double> >,hash<unsigned int>,equal_to<unsigned int> > PCFG_NTmap;
typedef google::dense_hash_map<string,set<pair<unsigned int, double> >,hash<string>,equal_to<string> > PCFG_PTmap;

struct PData {
    
    PData(unsigned int sym_) : sym(sym_) {
        in = 0;
        out = 0;
    }

    unsigned int sym;
    double in;
    double out;
    double prob;
    vector<pair<PData*, double> > ukids;
    vector<pair< pair<PData*, PData*>, double> > bkids;
    
};

typedef google::dense_hash_map<unsigned int,PData*,hash<unsigned int>,equal_to<unsigned int> > PCell;

struct PChart {

    unsigned int size;
    PCell** cells;
    
    PChart(unsigned int size_) : size(size_) {
        cells = new PCell*[size];
        for(size_t i=0;i<size;++i) {
            size_t rowS = size-i;
            cells[i] = new PCell[rowS];
            for(size_t j=0;j<rowS;++j) {
                cells[i][j].set_empty_key(100000);
                cells[i][j].set_deleted_key(100001);
            }
        }
    }

    ~PChart() {
        for(int i=0;i<size;++i) {
            size_t rowS = size-i;
            for(size_t j=0;j<rowS;++j) {
                PCell& cc = cells[i][j];
                for(PCell::iterator iter = cc.begin();iter != cc.end();++iter) {
                    delete iter->second;
                }
            }
            delete[] cells[i];
        }
        delete[] cells;
    }

    void print(string* syms) {
        for(int i=0;i<size;++i) {
            int rowS = size-i;
            for(int j=0;j<rowS;++j) {
                PCell& pc = cells[i][j];
                if(!pc.empty()) {
                    printf("CELL %d %d\n",i,j);
                    for(PCell::iterator iter = pc.begin();iter != pc.end();++iter) {
                        PData* d = iter->second;
                        printf("%s %E %E %E\n",syms[d->sym].c_str(),d->in,d->out,d->prob);
                    }
                }
            }
        }
    }

    vector<EvalItem> getEItems(string* syms) {
        vector<EvalItem> predict;
        for(int i=0;i<size;++i) {
            int rowS = size-i;
            for(int j=0;j<rowS;++j) {
                PCell& pc = cells[i][j];
                if(!pc.empty()) {
                    for(PCell::iterator iter = pc.begin();iter != pc.end();++iter) {
                        PData* d = iter->second;
                        EvalItem e(syms[d->sym],i,j,d->prob);
                        predict.push_back(e);
                    }
                }
            }
        }
        return predict;
    }
    
    void clean() {
        for(int i=0;i<size;++i) {
            int rowS = size-i;
            for(int j=0;j<rowS;++j) {
                PCell& pc = cells[i][j];
                if(!pc.empty()) {
                    for(PCell::iterator iter = pc.begin();iter != pc.end();++iter) {
                        PData* d = iter->second;
                        if(d->prob <= 0.0) {
                            pc.erase(iter->first);
                            delete d;
                        }
                    }
                }
            }
        }
    }
    
};

class PCFGParser {
public:
    PCFGParser(std::ifstream& ifs) {        
        string line;
        getline(ifs,line);
        nSym = strtol(line.c_str(),NULL,0);
        syms = new string[nSym];
        binaries = new PCFG_NTmap[nSym];
        unaries.set_empty_key(nSym+1);
        pterms.set_empty_key("EMPTYKEYSTRING");

        sym2base.set_empty_key("EMPTYKEY");

        for(int i=0;i<nSym;++i) {
            getline(ifs,line);
            string sym = line;
            syms[i] = sym;
            if(sym == "ROOT")
                rootSym = i;
            sym2base[sym] = i;
            binaries[i].set_empty_key(nSym+1);
        }

        getline(ifs,line);
        nRules = strtol(line.c_str(),NULL,0);
        printf("%d Rules\n",nRules);
        
        for(int i=0;i<nRules;++i) {
            getline(ifs,line);
            stringstream ss;
            ss << line;
            getline(ifs,line);
            double prob = atof(line.c_str());
            addRule(ss,prob,sym2base);
        }
    }

    ~PCFGParser() {
        if(syms != NULL) {
            delete[] syms;
            syms == NULL;
        }
        if(binaries != NULL) {
            delete[] binaries;
            binaries == NULL;
        }

    }

    void addRule(stringstream& ss, double p, S2Imap& sym2base) {
        char c;
        ss >> c;
        assert(c == '(');

        string sym;
        ss >> sym;
        assert(sym2base.count(sym) > 0);
        unsigned int lhs = sym2base[sym];
        pair<unsigned int,double> res = make_pair(lhs,p);
        vector<unsigned int> rhs;
        
        while(true) {
            c = ss.peek();
            if(c == ' ') { //if it's a space
                ss.ignore(1); //ignore that space
            } else if(c == '(') { //nonterminal
                ss >> c;
                ss >> sym;
                rhs.push_back(sym2base[sym]);
                while(c != ')') {
                    ss >> c;
                }
            } else if (c == ')') {
                if(rhs.size() == 1) {
                    PCFG_NTmap::iterator fter = unaries.find(rhs[0]);
                    if(fter == unaries.end()) {
                        set<pair<unsigned int, double> > myset;
                        myset.insert(res);
                        unaries[rhs[0]] = myset;
                    } else {
                        fter->second.insert(res);
                    }
                } else if (rhs.size() == 2) {
                    PCFG_NTmap& mymap = binaries[rhs[0]];
                    PCFG_NTmap::iterator fter = mymap.find(rhs[1]);
                    if(fter == mymap.end()) {
                        set<pair<unsigned int, double> > myset;
                        myset.insert(res);
                        mymap[rhs[1]] = myset;
                    } else {
                        fter->second.insert(res);
                    }
                } else {
                    printf("this is probably not a binarized pcfg - more than 2 RHS\n");
                    throw "bad rhs size";
                }
                
                return;
            } else { //terminal
                string term = "";
                while(ss.peek() != ')') {
                    term += ss.get();
                }

                //printf("%s\n",term.c_str());
                
                PCFG_PTmap::iterator fter = pterms.find(term);
                if(fter == pterms.end()) {
                    set<pair<unsigned int, double> > myset;
                    myset.insert(res);
                    pterms[term] = myset;
                } else {
                    fter->second.insert(res);
                }
                
                return;
            }
        }
    }

    void doUnaries(PChart* chart, size_t from, size_t to) {

        PCell& cell = chart->cells[from][to];
        
        for(unsigned int i=0;i<nSym;++i) { //assumes topological sort
            PCell::iterator fter = cell.find(i);
            if(fter != cell.end()) { //if the symbol is in the cell
                PData* ref = fter->second;
                PCFG_NTmap::iterator uter = unaries.find(i);
                if(uter != unaries.end()) {
                    for(set<pair<unsigned int, double> >::iterator iter = uter->second.begin();iter != uter->second.end();++iter) {
                        unsigned int sym = iter->first;
                        double p = iter->second;

                        PCell::iterator titer = cell.find(sym);
                        if(titer == cell.end()) { //make new data!
                            PData* d = new PData(sym);
                            d->in = ref->in*p;
                            d->ukids.push_back(make_pair(ref,p));
                            cell.insert(make_pair(sym,d));
                        } else { //add to previous data
                            PData* d = titer->second;
                            d->in += ref->in*p;
                            d->ukids.push_back(make_pair(ref,p));
                        }
                    }
                }
            }
        }

    }

    void doOutside(PData* d, double totalP) {

        for(size_t i=0;i<d->ukids.size();++i) {
            pair<PData*,double>& uk = d->ukids[i];
            PData* d2 = uk.first;
            //printf("UK : %s\n",syms[d2.sym].c_str());
            double p = uk.second;
            d2->out += d->out*p;
        }

        for(size_t i=0;i<d->bkids.size();++i) {
            pair<pair<PData*,PData*>,double>& uk = d->bkids[i];
            PData* d1 = uk.first.first;
            PData* d2 = uk.first.second;
            //printf("BK : %s %s\n",syms[d1.sym].c_str(),syms[d2.sym].c_str());
            double p = uk.second;
            d2->out += d->out*d1->in*p;
            d1->out += d->out*d2->in*p;
        }

        d->prob = d->out*d->in/totalP;
        
    }
    
    PChart* getChart(vector<string>& terms) {
        size_t nT = terms.size();
        
        PChart* chart = new PChart(terms.size());

        //insides
        for(size_t i=0;i<terms.size();++i) {

            PCFG_PTmap::iterator fter = pterms.find(terms[i]);
            if(fter == pterms.end()) {
                string s = (terms[i] + " not found");
                //printf("%s\n",s.c_str());
                throw s;
            } else {
                for(set<pair<unsigned int, double> >::iterator iter = fter->second.begin();iter != fter->second.end();++iter) {
                    unsigned int sym = iter->first;
                    PData* d = new PData(sym);
                    d->in = iter->second;

                    chart->cells[0][i].insert(make_pair(sym,d));
                }
                doUnaries(chart,0,i);
            }
        }

        //printf("CHART AFTER PTS\n");
        //chart->print(syms);
        
        for(size_t spanL=1;spanL<nT;++spanL) { //consider spans of length 2 to nT            
            for(size_t start=0;start<(nT-spanL);++start) { //starting at 0
                
                PCell& cc = chart->cells[spanL][start];
                
                for(size_t l1=0;l1<spanL;++l1) {
                    
                    size_t l2 = spanL - l1 - 1;
                    PCell& l = chart->cells[l1][start];
                    PCell& r = chart->cells[l2][start+l1+1];

                    if(l.empty() || r.empty())
                        continue;
                    
                    for(PCell::iterator mter = l.begin();mter != l.end();++mter) {
                        PData* lD = mter->second;
                        PCFG_NTmap& lmap = binaries[lD->sym];
                        //printf("looking at binary rules from %s\n",syms[lD.sym].c_str());
                        for(PCell::iterator nter = r.begin();nter != r.end();++nter) {
                            PData* rD = nter->second;
                            PCFG_NTmap::iterator fter = lmap.find(rD->sym);
                            if(fter != lmap.end()) {
                                //printf("found one with %s\n",syms[rD.sym].c_str());
                                for(set<pair<unsigned int, double> >::iterator iter = fter->second.begin();iter != fter->second.end();++iter) {
                                    unsigned int sym = iter->first;
                                    double p = iter->second;
                                    //printf("gives %s\n",syms[sym].c_str());
                                    PCell::iterator look = cc.find(sym);
                                    if(look == cc.end()) {
                                        PData* d = new PData(sym);
                                        d->in = p*lD->in*rD->in;
                                        d->bkids.push_back(make_pair(make_pair(lD,rD),p));
                                        cc[sym] = d;
                                    } else {
                                        PData* d = look->second;
                                        d->in += p*lD->in*rD->in;
                                        d->bkids.push_back(make_pair(make_pair(lD,rD),p));
                                    }
                                }
                            }
                        }
                    }
                }

                doUnaries(chart,spanL,start);
            }
        }

        //printf("\nCHART AFTER INSIDES\n");
        //chart->print(syms);
        //now the insides are full

        //find root sym

        PCell& cc = chart->cells[nT-1][0];
        PCell::iterator findo = cc.find(rootSym);
        if(findo == cc.end()) {
            //chart->print(syms);
            chart->clean();
            return NULL;
        }

        PData* d = findo->second;
        d->out = 1.0;
        double totalP = d->in;
        
        for(size_t i=0;i<nT;++i) {
            size_t len = nT-i-1;
            for(size_t start = 0;start < nT-len;++start) {
                PCell& ccell = chart->cells[len][start];
                for(int ss = nSym-1;ss>=0;--ss) {
                    unsigned int s = ss;
                    PCell::iterator findr = ccell.find(s);
                    if(findr != ccell.end()) {                    
                        PData* d = findr->second;
                        //printf("DO OUT - %s %d %d\n",syms[d.sym].c_str(),len,start);
                        doOutside(d,totalP);
                    }
                }
            }
        }

        chart->clean();
        
        //printf("\nCHART AT END\n");
        //chart->print(syms);
        
        return chart;
    }
    
    unsigned int nSym,rootSym,nRules;
    string* syms;
    S2Imap sym2base;
    PCFG_NTmap* binaries;
    PCFG_NTmap unaries;
    PCFG_PTmap pterms;

};

#endif
