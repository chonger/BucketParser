#ifndef BUCKET_EM
#define BUCKET_EM
#include "BucketGrammar.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <boost/math/special_functions/digamma.hpp>

using namespace std;

struct EMItem {

    EMItem(unsigned int sym_, vector<unsigned int> buckets_) :
        sym(sym_), buckets(buckets_), prob(0.0)
    {
        initH();
    }
    
    EMItem(unsigned int sym_, vector<unsigned int> buckets_,
           double prob_, vector<vector<EMItem*> >kids_) :
        sym(sym_), buckets(buckets_), prob(prob_), kids(kids_)
    {
        initH();
    }
    
    EMItem(unsigned int sym_,
           double prob_) :
        sym(sym_), buckets(), prob(prob_)
    {
        initH();
    }

        
    void initH() {
        hash<vector<unsigned int> > hf;
        hashv = sym ^ hf(buckets);
    }
    
    void add(double p,vector<EMItem*> k) {
        prob += p;        
        kids.push_back(k);
    }
    
    const unsigned int sym;
    const vector<unsigned int> buckets;
    double prob;
    vector<vector<EMItem*> > kids;
    size_t hashv;
    
private:
    EMItem(const EMItem &o) : sym(0), prob(0), hashv(0) {printf("!?!?!?!?!?!?\n");}
};

struct EMItemHash{
    size_t operator()(EMItem* const& k) const{
        return k->hashv;
    }
};

struct EMItemEq {
  bool operator()(EMItem* const& lhs, EMItem* const& rhs ) const
  {
      
      if(lhs->hashv == rhs->hashv) {
          equal_to<vector<unsigned int> > ef;
          return (lhs->sym == rhs->sym) && ef(lhs->buckets,rhs->buckets);
      } else
          return false;
      
  }
};

struct EMPEq {
  bool operator()(EMItem* const& lhs, EMItem* const& rhs ) const
  {
      return lhs == rhs;      
  }
};   

typedef google::dense_hash_set<EMItem*,EMItemHash,EMItemEq> EMap;
typedef google::dense_hash_set<EMItem*,EMItemHash,EMPEq> SeenMap;
typedef google::dense_hash_map<EMItem*,int,EMItemHash,EMPEq> EdgeMap;

struct EMLattice {

    EMLattice() {}
    
    EMLattice(unsigned int sym_, unsigned char nBucks_, vector<vector<unsigned int> > offs_) :
        inside(0.0),outside(0.0),sym(sym_),nBucks(nBucks_),offs(offs_){
        
    }    
    
    double inside,outside;
    unsigned int sym;
    unsigned char nBucks;
    vector<vector<unsigned int> > offs;
};


class BucketEM {
public:
    BucketEM(BucketGrammar& gr_, const char* trainF, unsigned int nBuckets_) : gr(gr_), nBuckets(nBuckets_) {

        gr.makePrior();
        
        sym2base.set_empty_key("EMPTYKEY");
        sym2base.set_deleted_key("DELETEDKEY");

        for(int i=0;i<gr.nSym;++i) {
            sym2base.insert(make_pair(gr.syms[i],i));
        }

        std::ifstream ifs(trainF);
        if(!ifs.is_open()) {
            printf("Invalid file at %s\n",trainF);
            exit(-2);
        }

        //gr.dump();
        
        //int c =0;
        while(ifs.good()){// && c++ < 1) {
            string treeS;
            getline(ifs,treeS);
            if(treeS.size() > 0)
                trainTreez.push_back(readTree(treeS));
        }
        ifs.close();
    }

    ~BucketEM() {
        for(vector<ParseTree*>::iterator iter = trainTreez.begin();iter != trainTreez.end();++iter) {
            delete *iter;
        }
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
            if(gr.tagP.find(l.sym) != gr.tagP.end())
                printf(" TAG ");
            if(gr.tsgP.find(l.sym) != gr.tsgP.end())
                printf(" TSG ");
            
            printf(" !!%f!!\n",kT);
        }
    }
    
    void recE(EMLattice& cur,vector<EMLattice>& lattice,
              double tot,I2Dmap& tsgE, I2Dmap& tagE, double* aE, double* nE) {
        

        //compute outsides and get expected counts
        
        double out = cur.outside;
        unsigned int s = cur.sym;

        if(cur.inside*out*tot - 1 > .1) {
            printLat(lattice);
            printf("myP %E \n",cur.inside*out*tot);
            throw "!";
        }
        
        if(cur.offs.size() == 0) {
            double myP = out*cur.inside;
            if(abs(myP*tot - 1) > .000001) {
                printLat(lattice);
                printf("MYP %f\n",myP);
                printf("TOT %f\n",1/tot);
                assert(abs(myP*tot-1) < .000001);
            }
        }
        
        for(vector<vector<unsigned int> >::iterator iter = cur.offs.begin();iter != cur.offs.end();++iter) {
            vector<unsigned int>& kids = *iter;

            assert(kids.size() > 0);
            assert(kids.size() <= 2);
            
            if(kids.size() == 1) {
                EMLattice& kid = lattice[kids[0]];
                unsigned int kS = kid.sym;
                if(cur.nBucks == kid.nBucks-1) { //going to a tag rule root
                    //involves deciding to adjoin and picking a tag rule
                    assert(cur.nBucks == kid.nBucks - 1);
                    //unsigned int bS = gr.baseSym[s];
                    double adjp = gr.adjoinProb(s);
                    double tag = gr.tagP[kS]; //TODO : Get this once!
                    double kOut = out*adjp*tag;
                    kid.outside += kOut;
                    double eCount = kOut * kid.inside * tot;
                    assert(eCount-1 < .1);
                    if(gr.aClassMap.count(s) > 0)
                        aE[gr.aClassMap[s]] += eCount;
                    tagE[kS] += eCount;
                } else {
                    if(s < gr.nSym) { //substitution - going to a tsg root node
                        //adjunction is disallowed here, so this only involves the tsg prob
                        double tsg = gr.tsgP[kS];
                        double kOut = out*tsg;
                        kid.outside += kOut;
                        double eCount = kOut * kid.inside * tot;
                        if(eCount - 1 > .1) {
                            printLat(lattice);
                            printf("EC %E \n",eCount);
                        }
                        assert(eCount-1 < .1);
                        tsgE[kS] += eCount;
                    } else {
                        if(s < gr.nSym*2) { //bc it's not a tsg root, it's a foot
                            //completely deterministic, no expectations to count
                            kid.outside += out;
                        } else { //normal unary
                            double adjp = gr.adjoinProb(s);
                            if(cur.nBucks < nBuckets && adjp > 0) {
                                //unsigned int bS = gr.baseSym[s];
                                double nadjp = 1-adjp;
                                double kOut = out*nadjp;
                                kid.outside += kOut;
                                double eCount = kOut * kid.inside * tot;
                                assert(eCount-1 < .1);
                                if(gr.aClassMap.count(s) > 0) {
                                    nE[gr.aClassMap[s]] += eCount;
                                }
                            } else {
                                kid.outside += out;
                            }
                        }
                    }
                }
            } else { //binary rule (must be in the transform)
                EMLattice& kidA = lattice[kids[0]];
                EMLattice& kidB = lattice[kids[1]];
                //                unsigned int bS = gr.baseSym[s];
                double adjp = gr.adjoinProb(s);
                if(cur.nBucks < nBuckets && adjp != 0) { //if this symbol isnt a tag root
                    //unsigned int bS = gr.baseSym[s];
                    double nadjp = 1-adjp;
                    double kOut = out*nadjp;
                    kidA.outside += kOut*kidB.inside;
                    kidB.outside += kOut*kidA.inside;
                    double eCount = kOut * kidA.inside * kidB.inside * tot;
                    assert(eCount-1 < .01);
                    if(gr.aClassMap.count(s) > 0)
                        nE[gr.aClassMap[s]] += eCount;
                } else {
                    kidA.outside += out*kidB.inside;
                    kidB.outside += out*kidA.inside;
                }
            }
        }
    }

    
    void recO1(EMItem* it, SeenMap& seen, EdgeMap& edges) {
        
        if(seen.count(it) > 0)
            return;
        seen.insert(it);

        for(vector<vector<EMItem*> >::iterator iter = it->kids.begin();iter != it->kids.end();++iter) {
            for(vector<EMItem*>::iterator jter = iter->begin();jter != iter->end();++jter) {
                edges[*jter] += 1;
                recO1(*jter,seen,edges);
            }
        }
        
    }

    void recO2(EMItem* it,vector<EMItem*>& ordered,EdgeMap& edges) {

        unsigned int count = edges[it];
        if(count <= 1) {
            ordered.push_back(it);
            for(vector<vector<EMItem*> >::iterator iter = it->kids.begin();iter != it->kids.end();++iter) {
                for(vector<EMItem*>::iterator jter = iter->begin();jter != iter->end();++jter) {
                    recO2(*jter,ordered,edges);
                }
            }
        } else {
            edges[it] = count-1;
        }
        
    }
    
    double emIter(bool vb) {

        I2Dmap tsgE,tagE;
        tsgE.set_empty_key(-1);
        tagE.set_empty_key(-1);
        double* aE = new double[gr.nAC];
        double* nE = new double[gr.nAC];
        for(int i=0;i<gr.nAC;++i) {
            aE[i] = 0.0;
            nE[i] = 0.0;
        }
        
        double ll = 0;                

        int c = 0;
        for(vector<ParseTree*>::iterator iter = trainTreez.begin();iter != trainTreez.end();++iter) {
            c++;
            if(c % 1000 == 0)
                printf(".\n");
            //printf("%s\n",gr.toString(*iter).c_str());
            
            EMItem* rootIter = doE(*iter);
            
            ll += log(rootIter->prob);
            //            printf("PROB = %e\n",rootIter->prob);
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

            //printf("%lu items in lattice\n",ordered.size());
            
            vector<EMLattice> lattice;

            EdgeMap inds;
            inds.set_empty_key(&eK);
            int x = 0;
            for(vector<EMItem*>::iterator j = ordered.begin();j != ordered.end();++j) {
                inds[*j] = x++;
            }
            
            
            for(size_t i = 0;i<ordered.size();++i) {

                EMItem* it = ordered[i];
                
                vector<vector< unsigned int> > offs;
                for(size_t j=0;j<it->kids.size();++j) {
                    vector<EMItem*>& cur = it->kids[j];
                    vector<unsigned int> newV;
                    for(size_t k=0;k<cur.size();++k) {
                        newV.push_back(inds[cur[k]]);
                    }
                    offs.push_back(newV);
                }
                EMLattice lat(it->sym,it->buckets.size(),offs);
                lat.inside = it->prob;
                lattice.push_back(lat);
            }
                        
            if(rootIter->prob <= 0) {
                printLat(lattice);
            }
            assert(rootIter->prob > 0);
            lattice[0].outside = 1.0;
            
            for(size_t i=0;i<lattice.size();++i) {
                //printf("!!!!!!!!!%lu\n",i);
                recE(lattice[i],lattice,tot,tsgE,tagE,aE,nE);
            }            

            for(vector<EMItem*>::iterator j = cleanme.begin();j != cleanme.end();++j) {
                delete *j;
            }
            cleanme.clear();

        }

        //M-STEP

        double adjPrior = .00000001;
        for(int i=0;i<gr.nAC;++i) {
            if(gr.adjoinP[i] > 0) {
                ll += adjPrior * log(gr.adjoinP[i] * (1-gr.adjoinP[i]));
            }
            if(aE[i]+nE[i] > 0 ) {                
                //                printf("NEW A %lu = %f/%f -> %f\n",i,aE[i],nE[i],(aE[i] + adjPrior) / (aE[i] + nE[i] + adjPrior*2));
                gr.adjoinP[i] = (aE[i] + adjPrior) / (aE[i] + nE[i] + adjPrior*2);
            }
        }

        double tsgTot[gr.nSym];
        double tsgPrior = 1;
        double tagTot[gr.nSym];
        double tagPrior = 1;
        for(size_t i=0;i<gr.nSym;++i) {
            tagTot[i] = 0.0;
            tsgTot[i] = 0.0;
        }


        //add counts from prior to expected counts
        for(I2Dmap::iterator iter = tagE.begin();iter!=tagE.end();++iter) {
            tagE[iter->first] = iter->second + gr.prior[iter->first] * tagPrior;
        }
        for(I2Dmap::iterator iter = tsgE.begin();iter!=tsgE.end();++iter) {
            tsgE[iter->first] = iter->second + gr.prior[iter->first] * tsgPrior;
        }

        //get sum of expected counts for each NT
        for(I2Dmap::iterator iter = gr.tagP.begin();iter!=gr.tagP.end();++iter) {
            unsigned int bS = gr.baseSym[iter->first];
            tagTot[bS] += tagE[iter->first];
        }
        for(I2Dmap::iterator iter = gr.tsgP.begin();iter!=gr.tsgP.end();++iter) {
            unsigned int bS = gr.baseSym[iter->first];
            double c = tsgE[iter->first];

            tsgTot[bS] += c;
        }

        if(vb) {
            for(size_t i=0;i<gr.nSym;++i) {
                if(tsgTot[i] > 0)
                    tsgTot[i] = boost::math::digamma(tsgTot[i]);
                if(tagTot[i] > 0)
                    tagTot[i] = boost::math::digamma(tagTot[i]);
            }        
            
            for(I2Dmap::iterator iter = gr.tsgP.begin();iter!=gr.tsgP.end();++iter) {
                unsigned int bS = gr.baseSym[iter->first];
                //if(iter->second > 0)
                //  ll += (gr.prior[iter->first] * tsgPrior-1) * log(iter->second);
                if(tsgTot[bS] > 0) {
                    
                    double newP = 0.0;
                    if(tsgE[iter->first] > 0)
                        newP = exp(boost::math::digamma(tsgE[iter->first]) - tsgTot[bS]);
                    /**
                       if(bS == 0) {
                       double p2 = tsgE[iter->first]/tsgTot[bS];
                       printf("P : V-%f N-%f\n",newP,p2);
                       }
                    */
                    //if(bS == 0) printf("NEWP %E - %E %E\n",newP,tsgE[iter->first],tsgTot[bS]);
                    gr.tsgP[iter->first] = newP;
                }
            }
        
            for(I2Dmap::iterator iter = gr.tagP.begin();iter!=gr.tagP.end();++iter) {
                //  if(iter->second > 0)
                //  ll += (gr.prior[iter->first] * tagPrior-1)* log(iter->second);
                unsigned int bS = gr.baseSym[iter->first];
                if(tagTot[bS] > 0) {
                    double c = tagE[iter->first];
                    if(c > 0)
                        gr.tagP[iter->first] = exp(boost::math::digamma(c) - tagTot[bS]);
                    else
                        gr.tagP[iter->first] = 0;
                }
            }
        } else {
            for(I2Dmap::iterator iter = gr.tsgP.begin();iter!=gr.tsgP.end();++iter) {
                unsigned int bS = gr.baseSym[iter->first];
                //if(iter->second > 0)
                //    ll += (gr.prior[iter->first] * tsgPrior-1) * log(iter->second);
                if(tsgTot[bS] > 0) {
                    double newP = tsgE[iter->first] / tsgTot[bS];
                    gr.tsgP[iter->first] = newP;
                }
            }
        
            for(I2Dmap::iterator iter = gr.tagP.begin();iter!=gr.tagP.end();++iter) {
                //  if(iter->second > 0)
                //  ll += (gr.prior[iter->first] * tagPrior-1)* log(iter->second);
                unsigned int bS = gr.baseSym[iter->first];
                if(tagTot[bS] > 0) 
                    gr.tagP[iter->first] = tagE[iter->first]/tagTot[bS];
            }
        }

        delete[] aE;
        delete[] nE;

        return ll;
    }
    
    EMItem* doE(ParseTree* tree) {
        unsigned int rootI = 0;
        vector<EMItem*> es = doE(tree->root,rootI,tree->terms);
        for(vector<EMItem*>::iterator iter = es.begin();iter != es.end();++iter) {
            if((*iter)->sym == gr.rootSym) 
                return *iter;
        }
        //FAIL
        return NULL;
    }

    vector<EMItem*> expandSameSyms(vector<EMItem*>& ems) {

        //printf("expand syms - %lu\n",ems.size());
        
        //do adjoin/complete/deref but not real unary
        
        EMap done;
        EMItem eK(-1,0.0);

        done.set_empty_key(&eK);
        //here the input will all be etree nodes - they can be completed and then adjoined.

        //all the inital nodes will be distinct - assertion checks this assumption
        for(vector<EMItem*>::iterator iter = ems.begin();iter != ems.end();++iter) {
            EMItem* em = *iter;
            assert(done.find(em) == done.end()); //make sure we dont need to combine things
            done.insert(em);
        }


        //first try to dereference tag rules
        //note that by disallowing root adjunction, a dereferenced node cannot be a tag root

        for(vector<EMItem*>::iterator iter = ems.begin();iter != ems.end();++iter) {
            EMItem* em = *iter;
            unsigned int s = em->sym;

            I2Dmap::iterator tagF = gr.tagP.find(s);
            if(tagF != gr.tagP.end()) {
                //printf("DEREF - %u\n",em->buckets.size());
                assert(em->buckets.size() > 0);
                assert(em->buckets.size() <= nBuckets);
                vector<EMItem*> v;
                v.push_back(em);
                vector<unsigned int> newB = em->buckets;
                unsigned int newS = newB.back();
                assert(gr.adjoinProb(newS) > 0);
                double p = tagF->second * gr.adjoinProb(newS) * em->prob;
                newB.pop_back();
                EMItem* newEM = new EMItem(newS,newB);
                EMap::iterator newF = done.find(newEM);
                if(newF == done.end()) {
                    newEM->add(p,v);
                    done.insert(newEM);
                    cleanme.push_back(newEM);
                } else {
                    delete newEM;
                    EMItem* oldEM = *newF;
                    oldEM->add(p,v);
                }
            } 
        }

        //now that we're dereferenced, we will never to deref
        //again need to b/c no root adjunction (to deref again we would need a tag root in a bucket)
        //
        //done now contains only original signed nodes, derefed nodes (original signed nodes that are not tag roots)
        //
        //so make another pass - try to tsg complete or adjoin - both lead to base Symbols that cannot be adjoined or completed
        //
        //
        //

        
        vector<EMItem*> derefed;
        for(EMap::iterator iter = done.begin();iter != done.end();++iter) {
            derefed.push_back(*iter);
        }
        
        for(vector<EMItem*>::iterator iter = derefed.begin();iter != derefed.end();++iter) {

            EMItem* em = *iter;
            unsigned int s = em->sym;
            
            double adjp = adjoinProb(em);
            
            I2Dmap::iterator tsgF = gr.tsgP.find(s);
            if(tsgF != gr.tsgP.end()) { //COMPLETE TSG
                double p = tsgF->second * em->prob;
                assert(em->buckets.size() == 0); //buckets must be empty
                vector<EMItem*> v;
                v.push_back(em);
                EMItem* newEM = new EMItem(gr.baseSym[s],em->buckets);
                EMap::iterator newF = done.find(newEM);
                if(newF == done.end()) {
                    newEM->add(p,v);
                    done.insert(newEM);
                    cleanme.push_back(newEM);
                } else {
                    delete newEM;
                    EMItem* oldEM = *newF;
                    oldEM->add(p,v);
                }
            }

            if(adjp > 0) { //ADJOIN
                double p = em->prob; //uses the same prob - this move is deterministic
                vector<EMItem*> v;
                v.push_back(em);
                vector<unsigned int> bucks = em->buckets;
                bucks.push_back(s);
                EMItem* newEM = new EMItem(gr.baseSym[s] + gr.nSym,bucks);
                EMap::iterator newF = done.find(newEM);
                if(newF == done.end()) {                    
                    newEM->add(p,v);
                    done.insert(newEM);
                    cleanme.push_back(newEM);
                } else {
                    EMItem* oldEM = *newF;
                    oldEM->add(p,v);
                    delete newEM;
                }
            }
        
        }

        vector<EMItem*> doneV;
        for(EMap::iterator iter = done.begin();iter != done.end();++iter) {
            EMItem* em = *iter;
            doneV.push_back(em);
        }
        return doneV;
    }
    
    double adjoinProb(EMItem* it) { 
        if(it->buckets.size() < nBuckets)
            return gr.adjoinProb(it->sym);
        else
            return 0.0;
    }
        
    vector<EMItem*> doE(TreeNode* node, unsigned int& tInd, vector<string>& terms) {
        vector<EMItem*> ret;
        assert(node->k.size() <= 2);
        if(node->k.size() == 0) {
            string term = terms[tInd];
            tInd++;
            S2Setmap::iterator pts = gr.preterms.find(term);
            assert(pts != gr.preterms.end()); //assert that the terminal is found
            for(set<unsigned int>::iterator iter = pts->second.begin();iter != pts->second.end();++iter) {
                unsigned int s = *iter;
                if(gr.baseSym[s] == node->sym) {
                    EMItem* it = new EMItem(s,1000.0);
                    cleanme.push_back(it);
                    ret.push_back(it);
                }
            }
        } else {
            if(node->k.size() == 1) { //unary
                vector<EMItem*> kE = doE(node->k[0],tInd,terms);
                for(vector<EMItem*>::iterator iter = kE.begin();iter != kE.end();++iter) {
                    EMItem* e = *iter;
                    Umap::iterator pars = gr.umap.find(e->sym);
                    if(pars == gr.umap.end()) {
                        //no completions...delete e?
                    } else {
                        //NOTE: this only hits once
                        for(set<unsigned int>::iterator jter = pars->second.begin();jter != pars->second.end();++jter) {
                            unsigned int us = *jter;
                            if(gr.baseSym[us] == node->sym) {
                                vector<EMItem*> r;
                                r.push_back(e);
                                EMItem* it = new EMItem(us,e->buckets);
                                double p = (1-adjoinProb(it))*e->prob;
                                it->add(p,r);
                                ret.push_back(it);
                                cleanme.push_back(it);
                            }
                        }
                    }
                }
            } else { //two kids
                vector<EMItem*> lE = doE(node->k[0],tInd,terms);
                vector<EMItem*> rE = doE(node->k[1],tInd,terms);
                for(vector<EMItem*>::iterator liter = lE.begin();liter != lE.end();++liter) {
                    EMItem* l = *liter;
                    for(vector<EMItem*>::iterator riter = rE.begin();riter != rE.end();++riter) {
                        EMItem* r = *riter;
                        Bmap::iterator pars = gr.bmap.find(make_pair(l->sym,r->sym));
                        if(pars == gr.bmap.end()) {
                            //no completions...delete e?
                        } else {
                            for(set<unsigned int>::iterator jter = pars->second.begin();jter != pars->second.end();++jter) {
                                unsigned int us = *jter;
                                if(gr.baseSym[us] == node->sym) {

                                    vector<EMItem*> rr;
                                    rr.push_back(l);
                                    rr.push_back(r);

                                    vector<unsigned int> bucks = l->buckets; //this is wasting computation!
                                    if(r->buckets.size() > 0)
                                        bucks = r->buckets;
                                    EMItem* it = new EMItem(us,bucks);
                                    double p = l->prob * r->prob * (1-adjoinProb(it));
                                    it->add(p,rr);
                                    cleanme.push_back(it);
                                    ret.push_back(it);
                                }
                            }
                        }
                    }
                }
            }
        }
        ret = expandSameSyms(ret);
        //printf("Done with %s - %lu items\n",gr.syms[gr.baseSym[node->sym]].c_str(),ret.size());
        assert(ret.size() > 0);
        return ret;
    }    
    
    ParseTree* readTree(string& s) {
        //printf("READING TREE %s\n",s.c_str());
        stringstream ss(s);
        return new ParseTree(readNode(ss),getTerms(s));
    }
    
    vector<string> getTerms(string& s) {
        vector<string> terms;
        //printf("EXAMINE %s\n",s.c_str());
        stringstream ss(s);
        string w;
        while(ss.good()) {
            ss >> w;
            if(w[0] != '(') {
                //printf("got %s\n",w.c_str());
                string tok = w.substr(0,w.find_first_of(')'));
                if(tok != "<>") {
                    //printf("--- %s\n",tok.c_str());
                    terms.push_back(tok);
                }
                
            }
        }
        return terms;
    }
            
    TreeNode* readNode(stringstream& ss) {
        char c;
        ss >> c;
        assert(c == '(');

        string sym = "";
        ss >> sym;
        unsigned int symI = sym2base[sym];
        
        vector<TreeNode*> kids;
        
        while(true) {
            c = ss.peek();
            if(c == ' ') { //it's a space
                ss.ignore(1); //ignore that space
            } else if(c == '(') { //nonterminal
                //get child node index
                TreeNode* kid = readNode(ss);
                kids.push_back(kid);                
            } else if (c == ')') {
                ss.ignore(1); //burn closing paren
                //      printf("Finished Rule with sym %s and %lu kids\n",sym.c_str(),kids.size()-1);

                
                while(kids.size() > 2) { //try to add glue rule
                    //  printf("K = %lu\n",kids.size());
                    TreeNode* r = kids.back();
                    kids.pop_back();
                    TreeNode* l = kids.back();
                    kids.pop_back();
                    unsigned int glueI = gr.nSym*2;
                    vector<TreeNode*> nKids;
                    nKids.push_back(l);
                    nKids.push_back(r);
                    
                    kids.push_back(new TreeNode(glueI,nKids));
                }

                return new TreeNode(symI,kids);
                
            } else { //terminal
                string term = "";
                while(ss.peek() != ')') {
                    term += ss.get();
                }
                //                printf("got terminal %s\n",term.c_str());
                ss.ignore(1); //burn closing paren
                vector<TreeNode*> noKids;
                return new TreeNode(symI,noKids);
            }
        }
        return NULL;
    }

    vector<EMItem*> cleanme;
    S2Imap sym2base;
    BucketGrammar& gr;
    vector<ParseTree*> trainTreez;
    unsigned int nBuckets;
};


#endif
