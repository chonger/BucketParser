#ifndef BUCKET_VB
#define BUCKET_VB

#include "VBGrammar.hpp"
#include "BucketEM.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <boost/math/special_functions/digamma.hpp>

typedef google::dense_hash_map<EMItem*,string,EMItemHash,EMPEq> TermMap;

class BucketVB {
public:
    
    BucketVB(VBGrammar& gr_) : gr(gr_) {
        
    }

    void printLat(vector<EMLattice>& lattice) {
        for(size_t i=0;i<lattice.size();++i) {
            EMLattice& l = lattice[i];
            double kT = 0;
            printf("%lu : %u(%s)/%u I=%e O=%e T=%e",i,l.sym,gr.syms[gr.baseSym[l.sym]].c_str(),l.nBucks,l.inside,l.outside,l.inside*l.outside);
            for(size_t j=0;j<l.offs.size();++j) {
                printf("[ ");
                vector<unsigned int>& o = l.offs[j];
                for(size_t k=0;k<o.size();++k){
                    printf("%u ",o[k]);
                    kT += lattice[o[k]].inside * lattice[o[k]].outside;
                }
                printf("]");
            }
            if(gr.tsgP.find(l.sym) != gr.tsgP.end())
                printf(" TSG ");
            
            printf(" !!%f!!\n",kT);
        }
    }

    void reformGrammar() {

        double cutoff = .001;
        double sampleCutoff = .001;
        
        T2Dmap treemap;
        treemap.set_empty_key(gr.eTree);

        double pLeft[gr.nSym];
        for(size_t i=0;i<gr.nSym;++i)
            pLeft[i] = 1.0;
        
        //eliminate rules that have a probability less than X (~10^-3 ?)
        for(T2Imap::iterator iter = gr.ruleTreez.begin();iter != gr.ruleTreez.end();++iter) {
            unsigned int ind = iter->second;
            unsigned int b = gr.baseSym[ind];
            double p = gr.tsgP[ind];
            if(p >= cutoff) {
                treemap[iter->first] = p;
                pLeft[b] -= p;
            } else {
                //printf("DROP : %s\n",iter->first->hashv.c_str());
            }
        }

        //DUMMY VARS FOR IO FUNCTION
        I2Dmap tsgE; //the expected counts of each tsg rule in the grammar
        tsgE.set_empty_key(-1);
        double alphaE[gr.nSym];
        for(size_t i=0;i<gr.nSym;++i) {
            alphaE[i] = 0.0;
        }
        
        //add rules by sampling from estimated trees
        for(vector<ParseTree*>::iterator iter = gr.trainTreez.begin();iter != gr.trainTreez.end();++iter) {
            ParseTree* tree = *iter;
            pair<vector<EMItem*>,double> res = insideOutside(tree,tsgE,alphaE);

            //tree->setHashString(gr.syms);
            // printf("SAMPLE %s\n",tree->hashv.c_str());
            
            double norm = 1.0 / res.first[0]->prob;
            
            for(size_t i=0;i<res.first.size();++i) {
                EMItem* cur = res.first[i];
                if(cur->sym < gr.nSym) {
                    //printf("Consider %s\n",gr.syms[cur->sym].c_str());
                    for(size_t j = 0;j<cur->kids.size();++j) {
                        vector<EMItem*>& kids = cur->kids[j].first;
                        if(kids.size() == 1 && kids[0]->sym == cur->sym + gr.nSym) { //its a base to a prime
                            double p = cur->oprob * cur->kids[j].second * norm;
                            if(p > ((double) rand()) / RAND_MAX && pLeft[cur->sym] > sampleCutoff) {
                                vector<string> newterms;
                                TreeNode* n = sampleTree(kids[0],newterms);
                                //       printf("NT\n");
                                //for(size_t k=0;k<newterms.size();++k) {
                                //    printf("%s\n",newterms[k].c_str());
                                //}
                                ParseTree* pt = new ParseTree(n,newterms);
                                pt->setHashString(gr.syms);
                                double stickW = .25; //TODO - sample this from a BETa
                                double myP = pLeft[cur->sym] * stickW;
                                treemap[pt] = myP;
                                pLeft[cur->sym] -= myP;
                                assert(pLeft > 0);
                            }    
                        }
                    }
                }
            }

            //cleanup
            for(vector<EMItem*>::iterator j = cleanme.begin();j != cleanme.end();++j) {
                delete *j;
            }
            cleanme.clear();
        }

        //how to initialize parameters for next EM round?
        vector<ParseTree*> tz;
        vector<double> pz;

        printf("New Grammar - %u rules\n",treemap.size());
        
        //        printf("NEW GRAMMAR\n");
        for(T2Dmap::iterator iter = treemap.begin();iter != treemap.end();++iter) {
            //   printf("%s %f\n",iter->first->hashv.c_str(),iter->second);
            tz.push_back(iter->first);
            pz.push_back(iter->second);
        }
        gr.makeTransform(tz,pz);
        /**
        printf("NEW ALPHAS\n");
        for(size_t i=0;i<gr.nSym;++i)
            printf("%s %f\n",gr.syms[i].c_str(),gr.alpha[i]);
        */
    }

    TreeNode* sampleTree(EMItem* cur, vector<string>& terms) {

        //        printf("Sample %s\n",gr.syms[cur->sym - gr.nSym].c_str());
        
        vector<TreeNode*> k;
        if(cur->kids.size() == 0) {
            string newT = *(cur->term);
            //   printf("ADD TERM %s\n",newT.c_str());
            terms.push_back(newT);
            return new TreeNode(cur->sym-gr.nSym,k);
        }
        
        //all kids are valid moves 
        //in a unary there will be 2, in a binary 4

        double r = ((double) rand()) / RAND_MAX;
        //printf("%f\n",r);
        double tot = 0.0;
        for(size_t i = 0;i<cur->kids.size();++i) {
            tot += cur->kids[i].second;
        }
        
        double sofar = 0.0;
        bool cont = true;
        for(size_t i = 0;i<cur->kids.size() && cont;++i) {
            double p = cur->kids[i].second / tot;
            sofar += p;
            if(sofar >= r) {
                cont = false;
                vector<EMItem*> choice = cur->kids[i].first;
                for(size_t j=0;j<choice.size();++j) {
                    EMItem* ch = choice[j];
                    if(ch->sym < gr.nSym) {
                        vector<TreeNode*> kk;
                        //            printf("ADD TERM <>\n");
                        terms.push_back("<>");
                        k.push_back(new TreeNode(ch->sym,kk));
                    } else {
                        k.push_back(sampleTree(ch,terms));
                    }
                }
                
            }
        }

        return new TreeNode(cur->sym-gr.nSym,k);
    }
    
    pair<vector<EMItem*>,double> insideOutside(ParseTree* pt,I2Dmap& tsgE,double* alphaE) {
        
        EMItem* rootIter = fillChart(pt);  //fill the chart and calculate inside probs
             
        double ll = log(rootIter->prob);
        //printf("PROB = %e\n",rootIter->prob);
        double tot = 1.0/rootIter->prob;
        
        //topological sort
        SeenMap seen;
        EdgeMap edges;
        EMItem eK(-1,0.0);
        seen.set_empty_key(&eK);
        edges.set_empty_key(&eK);
        vector<EMItem*> ordered;
        recO1(rootIter,seen,edges);
        recO2(rootIter,ordered,edges);

        ordered[0]->oprob = 1.0;
        for(size_t i = 0;i<ordered.size();++i) {
            EMItem* it = ordered[i]; 
            recE(it,tot,tsgE,alphaE);
        }

        return make_pair(ordered,ll);        
    }
    
    double emIter(bool vb) {

        I2Dmap tsgE; //the expected counts of each tsg rule in the grammar
        tsgE.set_empty_key(-1);
        double alphaE[gr.nSym];
        for(size_t i=0;i<gr.nSym;++i) {
            alphaE[i] = 1.0;
        }
        for(I2Dmap::iterator iter = gr.tsgP.begin();iter != gr.tsgP.end();++iter) {
            tsgE[iter->first] = 1.0;
        }
        
        double ll = 0; //data log likelihood
        for(vector<ParseTree*>::iterator iter = gr.trainTreez.begin();iter != gr.trainTreez.end();++iter) {
            pair<vector<EMItem*>,double> res = insideOutside(*iter,tsgE,alphaE);
            ll += res.second;
            //cleanup
            for(vector<EMItem*>::iterator j = cleanme.begin();j != cleanme.end();++j) {
                delete *j;
            }
            cleanme.clear();
        }

        //M-STEP
        double tsgTot[gr.nSym];
        for(size_t i=0;i<gr.nSym;++i) {
            //printf("ALPHA-E %s : %f\n",gr.syms[i].c_str(),alphaE[i]);
            tsgTot[i] = alphaE[i];
        }
        //get sum of expected counts for each NT
        for(I2Dmap::iterator iter = gr.tsgP.begin();iter!=gr.tsgP.end();++iter) {
            unsigned int bS = gr.baseSym[iter->first];
            double c = tsgE[iter->first];
            tsgTot[bS] += c;
        }

        for(I2Dmap::iterator iter = gr.tsgP.begin();iter!=gr.tsgP.end();++iter) {
            unsigned int bS = gr.baseSym[iter->first];
            double newP = tsgE[iter->first] / tsgTot[bS];
            //printf("TSG : %f\n",newP);
            gr.tsgP[iter->first] = newP;
        }
        for(size_t i=0;i<gr.nSym;++i) {
            gr.alpha[i] = alphaE[i]/tsgTot[i];
            //printf("ALPHA %s : %f\n",gr.syms[i].c_str(),gr.alpha[i]);
        }
        
        return ll;
    }
    

    
    //compute outsides and adds expected counts into I2Dmap
    void recE(EMItem* cur,
              double tot, //inverse of total probability
              I2Dmap& tsgE,//expected counts of TSG rules
              double* alphaE) 
    {
        double out = cur->oprob;
        unsigned int s = cur->sym;

        /**xs
        if(cur.offs.size() == 0) { //its a terminal
            //make sure inside * outside prob == tot^-1
            double myP = out*cur.inside;
            if(abs(myP*tot - 1) > .000001) {
                printLat(lattice);
                printf("MYP %f\n",myP);
                printf("TOT %f\n",1/tot);
                assert(abs(myP*tot-1) < .000001);
            }
        }
        */
        //for each set of possible children
        for(vector<pair<vector<EMItem*>,double> >::iterator iter = cur->kids.begin();iter != cur->kids.end();++iter) {

            vector<EMItem*>& kids = iter->first;

            assert(kids.size() > 0);
            assert(kids.size() <= 2);
            
            if(kids.size() == 1) { //unary rule
                EMItem* kid = kids[0];
                double p = out;
                unsigned int kS = kid->sym;
                if(gr.isPrime(s)) {
                    vector<unsigned int> key;
                    key.push_back(gr.baseSym[s]);
                    key.push_back(gr.baseSym[kS]);
                    if(gr.isPrime(kS)) {
                        p *= gr.pcfgNT[key] * (1-gr.cutProb[gr.baseSym[s]]); 
                    } else { //leaving prime
                        p *= gr.pcfgNT[key] * gr.cutProb[gr.baseSym[s]];
                    }
                } else {
                    if(gr.isPrime(kS)) {
                        p *= gr.alpha[gr.baseSym[s]];
                        //printf("%u !!!!!!! %f\n",s,p*kid.inside*tot);
                        alphaE[gr.baseSym[s]] += p * kid->prob * tot;
                    } else {
                        if(s < gr.nSym) { //substitution - going to a tsg root node
                            double tsgProb = gr.tsgP[kS];
                            //printf("TSG %f\n",tsgProb);
                            p *= tsgProb;

                            double eCount = p * kid->prob * tot;


                            assert(eCount-1 < .1);
                    
                            tsgE[kS] += eCount;
                        } else {
                            //tsg internal deterministic node, do nothing
                        }
                    }
                }
                kid->oprob += p;
            } else { //binary rule (must be in the transform)
                double p = out;
                EMItem* kidA = kids[0];
                EMItem* kidB = kids[1];
                unsigned int aSym = kidA->sym;
                unsigned int bSym = kidB->sym;
                if(gr.isPrime(s)) { 
                    vector<unsigned int> key;
                    key.push_back(gr.baseSym[s]);
                    key.push_back(gr.baseSym[aSym]);
                    key.push_back(gr.baseSym[bSym]);

                    if(aSym < gr.nSym) {
                        p *= gr.cutProb[aSym];
                    } else {
                        p *= 1.0 - gr.cutProb[gr.baseSym[aSym]];
                    }
                    if(bSym < gr.nSym) {
                        p *= gr.cutProb[bSym];
                    } else {
                        p *= 1.0 - gr.cutProb[gr.baseSym[bSym]];
                    }

                    //printf("!@!%f\n",gr.pcfgNT[key]);
                    p *= gr.pcfgNT[key]; 
                } else {
                    //deterministic
                }
                kidA->oprob += p*kidB->prob;
                kidB->oprob += p*kidA->prob;
            }
        }
    }

    void recO1(EMItem* it, SeenMap& seen, EdgeMap& edges) {
        
        if(seen.count(it) > 0)
            return;
        seen.insert(it);

        for(vector<pair<vector<EMItem*>,double > >::iterator iter = it->kids.begin();iter != it->kids.end();++iter) {
            for(vector<EMItem*>::iterator jter = iter->first.begin();jter != iter->first.end();++jter) {
                edges[*jter] += 1;
                recO1(*jter,seen,edges);
            }
        }
        
    }

    void recO2(EMItem* it,vector<EMItem*>& ordered,EdgeMap& edges) {

        unsigned int count = edges[it];
        if(count <= 1) {
            ordered.push_back(it);
            for(vector<pair<vector<EMItem*>,double > >::iterator iter = it->kids.begin();iter != it->kids.end();++iter) {
                for(vector<EMItem*>::iterator jter = iter->first.begin();jter != iter->first.end();++jter) {
                    recO2(*jter,ordered,edges);
                }
            }
        } else {
            edges[it] = count-1;
        }
        
    }
    
    
    void completeCell(unsigned int b, vector<EMItem*>& ems) {


        //printf("expand syms - %lu\n",ems.size());
        
        //1 - unprime the prime symbol
        //2 - complete TSGs


        EMItem* baseEM = new EMItem(b,0);
        cleanme.push_back(baseEM);
        
        for(vector<EMItem*>::iterator iter = ems.begin();iter != ems.end();++iter) {
            EMItem* em = *iter;
            unsigned int s = em->sym;
            if(gr.isPrime(s)) {
                //printf("??????%f\n",gr.alpha[gr.baseSym[s]]);
                double p = gr.alpha[gr.baseSym[s]] * em->prob;
                if(p <= 0) {
                    printf("%d\n",gr.baseSym[s]);
                    printf("A : %E\n",gr.alpha[gr.baseSym[s]]);
                    printf("B : %E\n",em->prob);
                }
                assert(p > 0);
                vector<EMItem*> v;
                v.push_back(em);
                baseEM->add(p,v);
            } else {
                I2Dmap::iterator tsgF = gr.tsgP.find(s);
                if(tsgF != gr.tsgP.end()) { //COMPLETE TSG
                    //printf("TSG IN %f\n",tsgF->second);
                    double p = tsgF->second * em->prob;
                    vector<EMItem*> v;
                    v.push_back(em);
                    baseEM->add(p,v);
                }
            }
        }

        ems.push_back(baseEM);

        //DIAGNOSIC
        for(size_t i=0;i<ems.size();++i) {
            if(ems[i]->prob <= 0) {
                printRec(ems[i]);
            }
        }
        
    }


    void printRec(EMItem* it) {
        string me = gr.syms[it->sym];
        double prob = it->prob;
        printf("%s %u %f\n",me.c_str(),it->sym,prob);
        for(size_t i=0;i<it->kids.size();++i) {
            printf("KIDS %d %f\n",i,it->kids[i].second);
            vector<EMItem*>& k = it->kids[i].first;
            for(size_t j=0;j<k.size();++j) {
                printf(" - %s %u\n",gr.syms[gr.baseSym[k[j]->sym]].c_str(),k[j]->sym);
            }
        }

    }

    //fills the chart and calculates inside probabilities
    EMItem* fillChart(ParseTree* tree) {
        unsigned int rootI = 0;
        vector<EMItem*> es = fillChart(tree->root,rootI,tree->terms);
        for(vector<EMItem*>::iterator iter = es.begin();iter != es.end();++iter) {
            if((*iter)->sym == gr.rootSym) 
                return *iter;
        }
        //FAIL
        return NULL;
    }

    /**
    void debug() {
        for(vector<EMItem*>::iterator kter = cleanme.begin();kter != cleanme.end();++kter) {
            EMItem* it = *kter;
            string ss = "";
        if(it->sym < gr.nSym) {
            ss =gr.syms[it->sym];
        }

        if(it->sym >= gr.nSym && it->sym < gr.nSym*2) {
            ss = gr.syms[it->sym-gr.nSym] + "*";
        }
        if(it->sym >= 2*gr.nSym)
            ss = gr.syms[gr.baseSym[it->sym]] + "'";
                
        
        printf("item %d - %s (%f)\n",it,ss.c_str(),it->prob);
            for(vector<vector<EMItem*> >::iterator iter = it->kids.begin();iter != it->kids.end();++iter) {
                for(vector<EMItem*>::iterator jter = iter->begin();jter != iter->end();++jter) {            
                    EMItem*  j = *jter;
                    printf("%d ",j);
                }
                printf("\n");
            }
        }
        int x;
        cin >> x;
    }
    */
    vector<EMItem*> fillChart(TreeNode* node, unsigned int& tInd, vector<string>& terms) {
        vector<EMItem*> ret;
        assert(node->k.size() <= 2);
        if(node->k.size() == 0) { //terminal
            string& term = terms[tInd++];
            pair<string,unsigned int> key = make_pair(term,node->sym);
            //            printf("term - %s\n",term.c_str());
            double primeProb = gr.pcfgPT[key] * (1.0 - gr.cutProb[node->sym]);
            EMItem* pit = new EMItem(gr.primeSym(node->sym),1000.0*primeProb);
            pit->term = &term;
            cleanme.push_back(pit);
            ret.push_back(pit);

            SIPmap::iterator pts = gr.preterms.find(key);
            if(pts != gr.preterms.end()) {
                unsigned int s = pts->second;
                EMItem* it = new EMItem(s,1000.0);
                //   printf("Found pterm %u\n",s);
                cleanme.push_back(it);
                ret.push_back(it);
            }

            
        } else {
            EMItem* primeEM = new EMItem(node->sym + gr.nSym,0);            
            
            if(node->k.size() == 1) { //unary
                
                vector<EMItem*> kE = fillChart(node->k[0],tInd,terms);

                for(vector<EMItem*>::iterator iter = kE.begin();iter != kE.end();++iter) {
                    EMItem* e = *iter;
                    unsigned int eSym = e->sym;
                    vector<unsigned int> key;
                    key.push_back(node->sym);
                    
                    if(gr.isPrime(eSym)) { //a prime rule must be used
                        key.push_back(gr.baseSym[eSym]);
                        double primeProb = gr.pcfgNT[key];
                        vector<EMItem*> r;
                        r.push_back(e);
                        double p = e->prob * primeProb * (1 - gr.cutProb[gr.baseSym[eSym]]);
                        assert(p > 0);
                        primeEM->add(p,r);
                    } else { //a TSG internal rule can be used
                        key.push_back(eSym);
                        V2Pmap::iterator pars = gr.rulemap.find(key);
                        if(pars != gr.rulemap.end()) {
                            //printf("TSG U FOUND\n");
                            unsigned int s = pars->second;
                            vector<EMItem*> r;
                            r.push_back(e);
                            EMItem* it = new EMItem(s,0);
                            double p = e->prob;
                            it->add(p,r);
                            ret.push_back(it);
                            cleanme.push_back(it);
                        }

                        if(eSym < gr.nSym) { //maybe do a prime rule
                            double p = e->prob * gr.pcfgNT[key] * gr.cutProb[gr.baseSym[eSym]];
                            vector<EMItem*> r;
                            r.push_back(e);
                            
                            if(p <= 0) {
                                printf("%s\n",gr.syms[eSym].c_str());
                                printf("fail to fill chart with unary prime\n");
                                printf("kid inside = %f\n",e->prob);
                                printf("pcfg prob = %f\n",gr.pcfgNT[key]);
                                printf("cut prob = %f\n",gr.cutProb[gr.baseSym[eSym]]);
                            }
                            assert(p > 0);
                            primeEM->add(p,r);
                        }
                    } 
                }
            } else { //two kids
                vector<EMItem*> lE = fillChart(node->k[0],tInd,terms);
                vector<EMItem*> rE = fillChart(node->k[1],tInd,terms);
                for(vector<EMItem*>::iterator liter = lE.begin();liter != lE.end();++liter) {
                    EMItem* l = *liter;
                    for(vector<EMItem*>::iterator riter = rE.begin();riter != rE.end();++riter) {
                        EMItem* r = *riter;
                        unsigned int rSym = r->sym;
                        unsigned int lSym = l->sym;

                        if(rSym < gr.nSym*2 && lSym < gr.nSym * 2) {
                            vector<unsigned int> key;
                            key.push_back(node->sym);
                            double p = l->prob * r->prob;
                            if(gr.isPrime(lSym)) {
                                key.push_back(lSym - gr.nSym);
                                p *= (1.0 - gr.cutProb[lSym-gr.nSym]);
                            } else {
                                key.push_back(lSym);
                                p *= gr.cutProb[lSym];
                            }
                            if(gr.isPrime(rSym)) {
                                key.push_back(rSym - gr.nSym);
                                p *= (1.0 - gr.cutProb[rSym-gr.nSym]);
                            } else {
                                key.push_back(rSym);
                                p *= gr.cutProb[rSym];
                            }
                            p *= gr.pcfgNT[key];
                            vector<EMItem*> rr;
                            rr.push_back(l);
                            rr.push_back(r);
                            primeEM->add(p,rr);
                        }

                        //check out the possibility of TSG internal rules
                        vector<unsigned int> key;
                        key.push_back(node->sym);
                        key.push_back(lSym);
                        key.push_back(rSym);
                        //printf("lookup binary rule %u -> %u %u\n",node->sym,lSym,rSym);
                        V2Pmap::iterator fter = gr.rulemap.find(key);
                        if(fter != gr.rulemap.end()) {
                            //printf("found binary rule %u(%u) -> %u %u\n",node->sym,fter->second,lSym,rSym);
                            vector<EMItem*> rr;
                            rr.push_back(l);
                            rr.push_back(r);
                            EMItem* it = new EMItem(fter->second,0);
                            double p = l->prob * r->prob;
                            it->add(p,rr);
                            cleanme.push_back(it);
                            ret.push_back(it);
                        }
                    }
                }
            }
            cleanme.push_back(primeEM);
            ret.push_back(primeEM);
        }
        completeCell(node->sym,ret);
        //printf("Done with %s - %lu items\n",gr.syms[gr.baseSym[node->sym]].c_str(),ret.size());
        assert(ret.size() > 0);
        return ret;
    }    
    
    vector<EMItem*> cleanme;
    VBGrammar& gr;

};


#endif
