#ifndef BUCKETPARSER
#define BUCKETPARSER 1

#include <cstddef>
#include <string>
#include <sstream>
#include <vector>
#include "BucketGrammar.hpp"
#include "Chart.hpp"
#include "PCFGParser.hpp"

using namespace std;

class BucketParser {
public:
    BucketParser(BucketGrammar& _gr, unsigned char nBuckets_) :
        gr(_gr), nBuckets(nBuckets_) {
        nullSym = new ChartItem(-1,0.0);
        delSym = new ChartItem(-2,0.0);
    }

    ~BucketParser() {
        if(nullSym != NULL)
            delete nullSym;
        nullSym = NULL;
        if(delSym != NULL)
            delete delSym;
        delSym = NULL;
    }

    string parse(string parseMe) {

        vector<string> terms;
        stringstream ss;
        ss << parseMe;
        while(ss.good()) {
            string t;
            ss >> t;
            terms.push_back(t);
        }
        return parse(terms);
        
    }

    ChartItem* tryInsert(size_t i, size_t j, ChartCell** chart, ChartItem* item) {
        bool inserted = false;
        ChartItem* ret = NULL;
            
        if(insertCheck(i,j,chart,item)) {   
            ChartCell& cc = chart[i][j];
            double iProb = item->prob;
            
            //printf("TRY INSERT %s\n",cistring(item).c_str());
            ChartCell::iterator fter = cc.find(item);
            
            if(fter == cc.end()) {
                cc.insert(item);
                ret = item;
                inserted = true;
                assert(ret->kids.size() > 0);
            } else {
                if(iProb > (*fter)->prob) {
                    ChartItem* prev = *fter;
                    prev->set(iProb,item->kids);
                    ret = prev;
                    assert(ret->kids.size() > 0);
                    delete item;
                    inserted = true;
                }
            }
        }
        
        if(!inserted)
            delete item;
        
        return ret;
    }

    void fillChart(ChartCell** chart, vector<string>& terms, bool viterbi) {

        prefill(terms);
        
        size_t nT = terms.size();
        //fill first row
        for(size_t i=0;i<nT;++i) {
            S2Setmap::iterator pts = gr.preterms.find(terms[i]);
            if(pts == gr.preterms.end()) {//assert that the terminal is found
                deleteChart(chart,terms);
                printf("Terminal %s not found\n",terms[i].c_str());
                return;// failString(terms);
            }
            for(set<unsigned int>::iterator iter = pts->second.begin();iter != pts->second.end();++iter) {
                ChartItem* item = new ChartItem(*iter,SCALE_FACTOR);
                chart[0][i].insert(item);
            }
            //printCell(0,i);
            doUnaries(0,i,chart);
            finalize(0,i,chart);
            
        }
        
        for(size_t spanL=1;spanL<nT;++spanL) { //consider spans of length 2 to nT

            //printf("Length %lu spans\n",spanL+1);
            
            for(size_t start=0;start<(nT-spanL);++start) { //starting at start

                ChartCell& cc = chart[spanL][start];
                
                for(size_t l1=0;l1<spanL;++l1) {
                    size_t l2 = spanL - l1 -1;

                    ChartCell& l = chart[l1][start];
                    ChartCell& r = chart[l2][start+l1+1];
                    vector<ChartItem*>& goodL = l.goodL;
                    vector<ChartItem*>& goodR = r.goodR;
                    Rmap& rlook = r.rlook;
                    size_t nR = rlook.size();
                    
                    for(vector<ChartItem*>::iterator lter = goodL.begin();lter != goodL.end();++lter) {
                        ChartItem* litem = *lter;
                        double naL = (1.0 - gr.adjoinProb(litem,nBuckets)) * litem->prob;
                        unsigned int lSym = litem->sym;
                        vector<pair<unsigned int, unsigned int> > llook = gr.leftlook[lSym]; //we know it has some L
                    
                        //see if there are fewer possible rules than RHSs
                        if(llook.size() < nR) { //if so, check every rule against RHSs
                            for(vector<pair<unsigned int, unsigned int> >::iterator rLter = llook.begin();rLter != llook.end();++rLter) {
                                unsigned int lC = rLter->first;
                                Rmap::iterator f = rlook.find(lC);
                                if(f != rlook.end()) { 
                                    for(vector<ChartItem*>::iterator rter = f->second.begin();rter != f->second.end();++rter) {
                                        ChartItem* ritem = *rter;
                                        double naR = (1.0 - gr.adjoinProb(ritem,nBuckets)) * ritem->prob;
                                        vector<unsigned int> buckets;
                                        if(litem->buckets.size() > 0)
                                            buckets = litem->buckets;
                                        if(ritem->buckets.size() > 0)
                                            buckets = ritem->buckets;
                                        unsigned int newSym = rLter->second;
                                        vector<ChartItem*> k;
                                        k.push_back(litem);
                                        k.push_back(ritem);
                                        double newProb = naL * naR;
                                        ChartItem* newItem = new ChartItem(newSym,buckets,newProb,k);
                                        tryInsert(spanL,start,chart,newItem);
                                    }
                                }
                            }
                        } else { //check possible RHS's in hashmap
                            for(vector<ChartItem*>::iterator rter = goodR.begin();rter != goodR.end();++rter) {
                                ChartItem* ritem = *rter;
                                Bmap::iterator bter = gr.bmap.find(make_pair(litem->sym,ritem->sym));
                                if(bter != gr.bmap.end()) {
                                    double naR = (1.0 - gr.adjoinProb(ritem,nBuckets)) * ritem->prob;                            
                                    for(set<unsigned int>::iterator iter = bter->second.begin();iter != bter->second.end();++iter) {
                                        vector<unsigned int> buckets;
                                        if(litem->buckets.size() > 0)
                                            buckets = litem->buckets;
                                        if(ritem->buckets.size() > 0)
                                            buckets = ritem->buckets;
                                        unsigned int newSym = *iter;
                                        vector<ChartItem*> k;
                                        k.push_back(litem);
                                        k.push_back(ritem);
                                        double newProb = naL * naR;
                                        ChartItem* newItem = new ChartItem(newSym,buckets,newProb,k);
                                        tryInsert(spanL,start,chart,newItem);
                                    }
                                }
                            }
                        }
                    }
                    
                }


                doUnaries(spanL,start,chart);
                
                finalize(spanL,start,chart);
                
                for(CCMap::iterator lter = cc.begin();lter != cc.end();++lter) {
                    assert((*lter)->kids.size() > 0);
                }

                //printCell(spanL,start);            
            }
        }

        postfill(terms);
    }

    ChartCell** makeChart(vector<string>& terms) {
        //build chart
        size_t nT = terms.size();
        ChartCell** chart = new ChartCell*[nT];
        for(size_t i=0;i<nT;++i) {
            chart[i] = new ChartCell[nT-i];
            for(size_t j=0;j<nT-i;++j) {
                chart[i][j].set_empty_key(nullSym);
                chart[i][j].set_deleted_key(delSym);
            }
        }
        return chart;
    }

    void deleteChart(ChartCell** chart, vector<string>& terms) {
        size_t nT = terms.size();
        size_t tot = 0;
        for(size_t i=0;i<nT;++i) {
            for(size_t j=0;j<nT-i;++j) {
                ChartCell& cc = chart[i][j];
                tot += cc.size();
                for(ChartCell::iterator iter = cc.begin();iter != cc.end();++iter) {
                    delete *iter;
                }
            }
            delete[] chart[i];
        }
        delete[] chart;
        printf("Contained %lu items\n",tot);
    }
    
    string parse(vector<string>& terms) {

        printf("Parsing, length %lu\n",terms.size());
        
        ChartCell** chart = makeChart(terms);
        fillChart(chart,terms,true);
        
        string parsed = "";
        size_t nT = terms.size();
        ChartItem root(gr.rootSym,1.0);
        ChartCell& rootCC = chart[nT-1][0];
        ChartCell::iterator fter = rootCC.find(&root);
        if(fter != rootCC.end()) {
            vector<string>::iterator t = terms.begin();
            parsed = getParse(*fter,t);
        } else {
            parsed = failString(terms);
            printf("\naFAIL %s\n",parsed.c_str());
        }
        
        deleteChart(chart,terms);
        
        return parsed;
    }

    string failString(vector<string>& terms) {
        stringstream ret;
        ret << "(ROOT";
        for(size_t i=0;i<terms.size();++i) {ret << " (X " + terms[i] + ")";}
        ret << ")";
        return ret.str();
    }

    string getParse(ChartItem* item, vector<string>::iterator& t) {

        unsigned int sym = item->sym;
        
        unsigned int bSym = gr.baseSym[item->sym];
        string symS = gr.syms[bSym];
        //printf("%s\n",symS.c_str());
        bool isB = symS.find('@') == 0 || symS.find('+') == 0 || symS == "GLUE";
        //bool isB = false;
        if(sym == bSym) {
            assert(item->kids.size() == 1);
            return getParse(item->kids[0],t);
        }
        if(gr.tagP.count(sym) > 0 || isB) {

            stringstream ret;
        
            ret << getParse(item->kids[0],t);

            for(size_t i=1;i<item->kids.size();++i) {
                ret << " " << getParse(item->kids[i],t);
            }

            //ret << ")";
        
            return ret.str();
            
        } else {
            
            stringstream ret;
        
            ret << "(" << gr.syms[bSym];

            if(item->kids.size() > 0) {            
                for(size_t i=0;i<item->kids.size();++i) {
                    ret << " " << getParse(item->kids[i],t);
                }
            } else {
                ret << " " + *t;
                ++t;
            }
            ret << ")";
        
            return ret.str();
        }
    }
    /**    
    void printCell(size_t i,size_t j) {
        printf("\nCELL %lu %lu\n",i,j);
        printCell(chart[i][j]);
    }

    void printCell(ChartCell& cc) {
        
        for(ChartCell::iterator iter = cc.begin();iter != cc.end();++iter) {
            ChartItem* it = *iter;

            vector<string> c;
            c.push_back("word");
            vector<string>::iterator ct = c.begin();
            printf("%s %e %s\n",cistring(it).c_str(),it->prob,getParse(it,ct).c_str());
   
            printf("%s %e\n",cistring(it).c_str(),it->prob);
        }
        
    }
**/
    
    string cistring(ChartItem* it) {
        stringstream stackS;
        stackS << it->sym  << " " << gr.syms[gr.baseSym[it->sym]] << " ";

        vector<unsigned int> s = it->buckets;
        //printf("NB : %u\n",it.buckets.size());
        while(!s.empty()) {
            unsigned int sym = s.back();
            s.pop_back();
            stackS << "(" << sym << " " << gr.syms[gr.baseSym[sym]] << ") ";
        }

        return stackS.str();
    }

    
    void doUnaries(size_t i, size_t j, ChartCell** chart) {

        ChartCell& cc = chart[i][j];
        
        vector<ChartItem*> proc;
        for(ChartCell::iterator iter = cc.begin();iter != cc.end();++iter) {
            proc.push_back(*iter);
        }
        
        while(!proc.empty()) {

            ChartItem* it = proc.back();
            double prob = it->prob;
            proc.pop_back();

            unsigned int sym = it->sym;

            //4 possibilities - adjunction, completion (TAG/TSG) , realUnary

            double adjP = gr.adjoinProb(it,nBuckets);
            double nadjP = 1-adjP;

            //COMPLETION
            //for tsg completion -> base sym, for tag completion -> top bucket
            I2Dmap::iterator tagF = gr.tagP.find(sym);
            if(tagF != gr.tagP.end()) { //TSG completion

                vector<unsigned int> buckets = it->buckets;
                unsigned int newSym = buckets.back();
                buckets.pop_back();
                vector<ChartItem*> k;
                k.push_back(it);
                double p = tagF->second;
                double newProb = prob*(nadjP)*p;
                ChartItem* newIt = new ChartItem(newSym,buckets,newProb,k);
                newIt = tryInsert(i,j,chart,newIt);
                if(newIt != NULL) {
                    proc.push_back(newIt);
                }
                
            } else {
                I2Dmap::iterator tsgF = gr.tsgP.find(sym);
                if(tsgF != gr.tsgP.end()) { //TAG completion
                    unsigned int newSym = gr.baseSym[sym];
                    //no change to buckets
                    vector<ChartItem*> k;
                    k.push_back(it);
                    double p = tsgF->second;
                    double newProb = prob*(nadjP)*p;
                    ChartItem* newIt = new ChartItem(newSym,it->buckets,newProb,k);
                    newIt = tryInsert(i,j,chart,newIt);
                    if(newIt != NULL) {
                        proc.push_back(newIt);
                    }

                } 

                //ADJUNCTION
                if(adjP > 0.0) {
                    unsigned int adjSym = gr.baseSym[sym] + gr.nSym;
                    vector<unsigned int> buckets = it->buckets;
                    buckets.push_back(sym);
                    
                    vector<ChartItem*> k;
                    k.push_back(it);
                    double newProb = prob*adjP;
                    ChartItem* newIt = new ChartItem(adjSym,buckets,newProb,k);
                    newIt = tryInsert(i,j,chart,newIt);
                    if(newIt != NULL) {
                        proc.push_back(newIt);
                    }
                }
            }

            //REAL UNARY
            
            Umap::iterator fter = gr.umap.find(sym);
            
            if(fter != gr.umap.end()) {
                for(set<unsigned int>::iterator iter = fter->second.begin();iter != fter->second.end();++iter) {
                    vector<unsigned int> buckets = it->buckets;
                    unsigned int newSym = *iter;
                    vector<ChartItem*> k;
                    k.push_back(it);
                    double newProb = prob*(nadjP);
                    ChartItem* newIt = new ChartItem(newSym,buckets,newProb,k);
                    newIt = tryInsert(i,j,chart,newIt);
                    if(newIt != NULL) {
                        proc.push_back(newIt);
                    }
                }
            }
        }

        //printCell(cc);
        
    }

    virtual bool insertCheck(size_t i, size_t j, ChartCell** chart, ChartItem* it) {
        return true;
    }
    
    virtual void finalize(size_t i, size_t j, ChartCell** chart) {
        ChartCell& cc = chart[i][j];
        double best = 0;
        for(CCMap::iterator lter = cc.begin();lter != cc.end();++lter) {
            double p = (*lter)->prob;
            if(p > best)
                best = p;
        }
        double cut = best * BEAM_WIDTH;
        
        for(CCMap::iterator lter = cc.begin();lter != cc.end();++lter) {
            unsigned int s = (*lter)->sym;
            double p = (*lter)->prob;
            if(p >= cut) {
                if(gr.canL[s])
                   cc.goodL.push_back(*lter);
                if(gr.canR[s]) {
                    cc.rlook[s].push_back(*lter);
                    cc.goodR.push_back(*lter);
                }
            } else {
                ChartItem* it = *lter;
                cc.erase(it);
                delete it;
            }
        }

        cc.resize(0);
    }

    virtual void prefill(vector<string>& terms) {
        return;
    }

    virtual void postfill(vector<string>& terms) {
        return;
    }

protected:
    BucketGrammar& gr; //the primary grammar
    
private:

    const static double SCALE_FACTOR = 1000.0;
    const static double BEAM_WIDTH = .00000000001;//.000000000000000000001;

    unsigned char nBuckets;
                   
    ChartItem* nullSym;
    ChartItem* delSym;
    
};



class CTFBucketParser : public BucketParser {
public:
    CTFBucketParser(BucketGrammar& gr, unsigned char nBuckets, PCFGParser* pp) : BucketParser(gr,nBuckets), pcfgParser(pp) {        
        
    }

    void prefill(vector<string>& terms) {
        
        ctfchart = pcfgParser->getChart(terms);

    }

    void postfill(vector<string>& terms) {

        delete ctfchart;
        ctfchart = NULL;
        
    }

    bool insertCheck(size_t i, size_t j, ChartCell** chart, ChartItem* it) {

        double cut = .0001;
        
        PCell& cc = ctfchart->cells[i][j];

        unsigned int s = it->sym;
        unsigned int bS = gr.baseSym[s];
        if(bS == gr.nSym*2)
            return true;
        if(bS > gr.nSym)
            bS -= gr.nSym;
        string baseS = gr.syms[bS];
        unsigned int ctfID = pcfgParser->sym2base[baseS];

        PCell::iterator findo = cc.find(ctfID);
        if(findo != cc.end() && findo->second->prob >= cut) {
            return true;
        }
        return false;
    }
    
private:
    PChart* ctfchart;
    PCFGParser* pcfgParser; 
};

#endif
