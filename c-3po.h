// [min,max]-3PO algorithm for constructing 2-hop cover
#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#define DEBUG(X) //cout << X << endl

// toggles for algorithmic features
#define MIN_MAX_CENTRALITY
#define UPDATE_CENTRALITY
#define RANDOMIZE
#define ESTIMATE_ANC_DESC

typedef uint32_t NodeID;

static const uint32_t DEFAULT = UINT32_MAX;
static const uint32_t DELETED = UINT32_MAX - 1;

struct Neighbor
{
    NodeID node;
#ifdef ESTIMATE_ANC_DESC
    uint32_t estimate;
    bool operator<(Neighbor other);
#endif
    Neighbor(NodeID node);
};

// vector of nodes with lazy deletion
class LazyList: public std::vector<Neighbor>
{
    size_type lazy_size = 0; // track actual data size
public:
    size_type size() const;
    void push_back(NodeID val);
    // notify of size reduction
    void operator--();
    // remove deleted nodes while preserving node order (call before iterating over nodes => no complexity increase)
    void shrink(const std::vector<uint32_t> &status, size_t start = 0);
    void clear();
    // erase non-deleted nodes given as sorted vector
    size_t erase_all(const std::vector<NodeID> &sorted);
};

// one-directional graph
class PartialGraph
{
public:
    std::vector<LazyList> neighbors;
    std::vector<NodeID> *last_visited;

    PartialGraph(size_t size, std::vector<NodeID> *last_visited);
    PartialGraph(const PartialGraph &g);
    bool deleted(NodeID node) const;
};

class DiGraph
{
    // for tracking deleted & visited nodes during propagation
    std::vector<NodeID> last_visited;
public:
    PartialGraph forward;
    PartialGraph backward;
#ifdef ESTIMATE_ANC_DESC
    std::vector<uint32_t> anc_estimate; // total #ancestors
    std::vector<uint32_t> desc_estimate; // total #descendants
#endif

    DiGraph(size_t size);
    DiGraph(const DiGraph &g);
    void insert_edge(NodeID from, NodeID to);
    void remove_node(NodeID node);
    uint64_t centrality(NodeID node) const;
    size_t size() const;
#ifdef ESTIMATE_ANC_DESC
    void estimate_anc_desc();
#endif
};

typedef std::vector<NodeID> LabelSet;

class TwoHopCover
{
public:
    std::vector<LabelSet> in;
    std::vector<LabelSet> out;

    TwoHopCover(size_t nodes);
    size_t labels() const;
    size_t self_labels() const;
    void map(std::vector<NodeID> f);
    // functions for reachability testing
    void erase_self();
    void insert_self();
    bool can_reach_all_self(NodeID from, NodeID to) const;
    bool can_reach_remote(NodeID from, NodeID to) const;
    bool can_reach_no_self(NodeID from, NodeID to) const;
};

TwoHopCover pick_propagate_prune(DiGraph &g, std::vector<NodeID> &pick_order);
DiGraph read_graph(std::istream &in);
DiGraph read_graph_from_file(char* filename);

// -------------------------- output functions --------------------------------

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> &v)
{
    os << "[";
    for ( size_t i = 0; i < v.size(); i++ )
    {
        os << (i ? ", " : " ") << i << ": " << v[i];
    }
    return os << " ]";
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2> &p)
{
    return os << "(" << p.first << ", " << p.second << ")";
}

std::ostream& operator<<(std::ostream &os, const Neighbor &n);
std::ostream& operator<<(std::ostream &os, const LazyList &l);
std::ostream& operator<<(std::ostream &os, const DiGraph &g);
std::ostream& operator<<(std::ostream &os, const TwoHopCover &c);
