// Implementation of [min,max]-3PO algorithm for constructing 2-hop cover
// compile with g++ -o c3po -O3 -Wall c-3po.cpp

#include <algorithm>
#include <cstdint>
#include <queue>
#include <vector>

using namespace std;

typedef uint32_t NodeID;

// vector of nodes with lazy deletion
class LazyList: public vector<NodeID>
{
    size_type lazy_size = 0; // track actual data size
public:

    size_type size() const // override (non-virtual)
    {
        return lazy_size;
    }

    void push_back(const NodeID &val) // override (non-virtual)
    {
        vector::push_back(val);
        ++lazy_size;
    }

    // notify of size reduction
    void operator--()
    {
        --lazy_size;
    }

    // remove deleted nodes
    void prune(const vector<bool> &deleted)
    {
        if ( lazy_size < vector::size() )
        {
            auto it = begin();
            while ( it != end() )
            {
                if ( deleted[*it] )
                {
                    *it = back();
                    pop_back();
                }
                else
                    ++it;
            }
        }
    }

    void clear() // override (non-virtual)
    {
        vector::clear();
        lazy_size = 0;
    }
};

// one-directional graph
class PartialGraph
{
public:
    vector<LazyList> neighbors;
    const vector<bool> &deleted;

    PartialGraph(const vector<bool> &deleted, size_t size) : neighbors(size), deleted(deleted) {}
};

class DiGraph
{
public:
    PartialGraph forward;
    PartialGraph backward;
    vector<bool> deleted;

    DiGraph(size_t size) : forward(deleted, size), backward(deleted, size), deleted(size) {}

    void remove_node(NodeID node)
    {
        for ( NodeID neighbor : forward.neighbors[node] )
            --backward.neighbors[neighbor];
        for ( NodeID neighbor : backward.neighbors[node] )
            --forward.neighbors[neighbor];
        forward.neighbors[node].clear();
        backward.neighbors[node].clear();
        deleted[node] = true;
    }

    uint64_t centrality(NodeID node)
    {
        uint64_t in = backward.neighbors[node].size();
        uint64_t out = forward.neighbors[node].size();
        return in <= out ? (in << 32) & out : (out << 32) & in;
    }

    size_t size()
    {
        return deleted.size();
    }
};

// break centrality-ties using pseudo-random permutation of node ID
class WeightedNode
{
public:
    uint64_t weight;
    uint32_t encoded;

    WeightedNode(uint64_t weight, uint32_t node) : weight(weight)
    {
        static const uint32_t mult = 2971215073u; // largest 32-bit Fibbonacci number
        encoded = node * mult;
    }

    uint32_t node()
    {
        static const uint32_t inverse = 1051402017u; // multiplicative inverse modulo 2^32
        return encoded * inverse;
    }

    bool operator<(WeightedNode other) const
    {
        return weight < other.weight || (weight == other.weight && encoded < other.encoded);
    }
};

// labels are stored as sorted vector - fast for small sets
class LabelSet : public vector<NodeID>
{
public:
    LabelSet(size_t capacity) : vector(capacity) {}

    bool contains(NodeID node) const
    {
        return std::binary_search(cbegin(), cend(), node);
    }

    void insert(NodeID node)
    {
        vector::insert(std::upper_bound(cbegin(), cend(), node), node);
    }

    bool intersects(const LabelSet &other) const
    {
        auto a = cbegin();
        auto b = other.cbegin();
        while ( a < cend() && b < other.cend() )
        {
            if ( *a < *b )
                ++a;
            else if ( *a > *b )
                ++b;
            else
                return true;
        }
        return false;
    }
};

struct TwoHopCover
{
    vector<LabelSet> in;
    vector<LabelSet> out;
    TwoHopCover(size_t nodes, size_t label_capacity=4) : in(nodes, label_capacity), out(nodes, label_capacity) {}
};


void propagate_prune(PartialGraph &g, NodeID node, LabelSet &node_labels, vector<LabelSet> &labels)
{
    // insert self label
    if ( g.neighbors[node].size() )
        node_labels.insert(node);
    // propagate remote labels
    vector<NodeID> stack;
    for ( NodeID neighbor : g.neighbors[node] )
        // will be removed after, so just check instead of pruning
        if ( !g.deleted[neighbor] )
            stack.push_back(neighbor);
    while ( !stack.empty() )
    {
        NodeID next = stack.back(); stack.pop_back();
        if ( !labels[next].intersects(labels[node]) )
        {
            labels[next].insert(node);
            // prune while we're here
            g.neighbors[next].prune(g.deleted);
            for ( NodeID neighbor : g.neighbors[next] )
                stack.push_back(neighbor);
        }
    }
}

TwoHopCover pick_propagate_prune(DiGraph &g)
{
    TwoHopCover labels(g.size());
    // order nodes by centrality
    priority_queue<WeightedNode> q;
    for ( NodeID node = 0; node < g.size(); ++node )
        q.push(WeightedNode(g.centrality(node), node));
    // pick top-ranked node, updating centrality lazily
    while ( !q.empty() )
    {
        WeightedNode weighted = q.top(); q.pop();
        NodeID node = weighted.node();
        uint64_t new_weight = g.centrality(node);
        if ( new_weight < weighted.weight )
        {
            weighted.weight = new_weight;
            q.push(weighted);
        }
        else
        {
            propagate_prune(g.forward, node, labels.out[node], labels.in);
            propagate_prune(g.backward, node, labels.in[node], labels.out);
            g.remove_node(node);
        }
    }
    return labels;
}

int main (int argc, char *argv[])
{
    return 0;
}
