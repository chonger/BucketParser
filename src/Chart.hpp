#ifndef BUCKETCHART
#define BUCKETCHART 1

#include <vector>
#include <google/sparse_hash_set>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <google/dense_hash_map>
#include <boost/functional/hash.hpp>
using namespace std;
using boost::hash;


struct ChartItem {
    
    ChartItem(unsigned int sym_, vector<unsigned int> buckets_,
              double prob_, vector<ChartItem*> kids_) :
        sym(sym_), buckets(buckets_), prob(prob_), kids(kids_)
    {
        //normal full on constuctor
        hash<vector<unsigned int> > hf;
        hashv = sym ^ hf(buckets);
    }

    ChartItem(unsigned int sym_,
              double prob_) :
        sym(sym_), prob(prob_)
    {
        hash<vector<unsigned int> > hf;
        hashv = sym ^ hf(buckets);
        //used for constucting preterminal chart items
    }

    void set(double p,vector<ChartItem*>&k) {
        prob = p;
        kids.clear();
        for(vector<ChartItem*>::iterator iter = k.begin();iter !=k.end();++iter) {
            kids.push_back(*iter);
        }
        assert(kids.size() > 0);
    }
    
    const unsigned int sym;
    const vector<unsigned int> buckets;
    double prob;
    vector<ChartItem*> kids;
    size_t hashv;
    
private:
    ChartItem(const ChartItem &o) : sym(0), prob(0), hashv(0) {}
};

struct ChartItemHash{
    size_t operator()(ChartItem* const& k) const{
        return k->hashv;
    }
};

struct ChartItemEq {
  bool operator()(ChartItem* const& lhs, ChartItem* const& rhs ) const
  {
      
      if(lhs->hashv == rhs->hashv) {
          equal_to<vector<unsigned int> > ef;
          return (lhs->sym == rhs->sym) && ef(lhs->buckets,rhs->buckets);
      } else
          return false;
      
  }
};   

typedef google::dense_hash_map<unsigned int,vector<ChartItem*> > Rmap;
typedef google::dense_hash_set<ChartItem*,ChartItemHash,ChartItemEq> CCMap;

class ChartCell : public CCMap {
public:

    ChartCell() {
        rlook.set_empty_key(-1);
    }
    
    Rmap rlook;
    vector<ChartItem*> goodL,goodR;

private:
    ChartCell(const ChartCell& o) {}
};



#endif
