// Implementation of [min,max]-3PO algorithm for constructing 2-hop cover
// compile with g++ -o c3po -O3 -Wall c-3po.cpp

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <queue>
#include <vector>
#include <sstream>
#include <string>

#define DEBUG(X) //cout << X << endl

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

    void insert_edge(NodeID from, NodeID to)
    {
        forward.neighbors[from].push_back(to);
        backward.neighbors[to].push_back(from);
        DEBUG("inserted edge " << from << " -> " << to);
    }

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
        return in <= out ? (in << 32) | out : (out << 32) | in;
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

class TwoHopCover
{
public:
    vector<LabelSet> in;
    vector<LabelSet> out;
    TwoHopCover(size_t nodes) : in(nodes), out(nodes) {}

    size_t labels() const
    {
        size_t sum = 0;
        for ( const LabelSet &ls : in )
            sum += ls.size();
        for ( const LabelSet &ls : out )
            sum += ls.size();
        return sum;
    }

    size_t self_labels() const
    {
        size_t sum = 0;
        for ( NodeID node = 0; node < in.size(); ++node )
        {
            if ( in[node].contains(node) )
                ++sum;
            if ( out[node].contains(node) )
                ++sum;
        }
        return sum;
    }

    bool can_reach(NodeID from, NodeID to) const
    {
        return from == to || out[from].intersects(in[to]);
    }
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
            DEBUG("picked " << node << " with C=[" << (new_weight >> 32) << ", " << (new_weight & 4294967295u) << "]");
            propagate_prune(g.forward, node, labels.out[node], labels.in);
            propagate_prune(g.backward, node, labels.in[node], labels.out);
            g.remove_node(node);
        }
    }
    return labels;
}

DiGraph read_graph()
{
    size_t nodes;
    string line;
    NodeID node;
    // first line contains #nodes (optionally followed by comment)
    cin >> nodes;
    DiGraph g(nodes);
    getline(cin, line);
    // remaining lines have the format node child_1 ... child_n
    while ( cin ) {
        getline(cin, line);
        // parse to integers
        istringstream ss(line);
        vector<NodeID> node_and_children;
        while ( ss >> node )
            node_and_children.push_back(node);
        // construct graph
        for ( size_t child = 1; child < node_and_children.size(); ++child )
            g.insert_edge(node_and_children[0], node_and_children[child]);
    }
    return g;
}

ostream& operator<<(ostream &os, const DiGraph &g);
ostream& operator<<(ostream &os, const TwoHopCover &c);

int main (int argc, char *argv[])
{
    DiGraph g = read_graph();
    DEBUG("g=" << g);
    auto t_start = chrono::high_resolution_clock::now();
    TwoHopCover cover = pick_propagate_prune(g);
    auto t_stop = chrono::high_resolution_clock::now();
    long dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
    DEBUG("cover=" << cover);
    // report stats
    size_t labels = cover.labels();
    size_t self_labels = cover.self_labels();
    cout << "found 2-hop cover with " << labels - self_labels << " remote and " << self_labels << " self labels in " << static_cast<double>(dur_ms) / 1000.0 << "s" << endl;
    return 0;
}

// -------------------------- output tools --------------------------

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

ostream& operator<<(ostream &os, const LazyList &l)
{
    return os << l.size() << " of " << static_cast<vector<NodeID>>(l);
}

ostream& operator<<(ostream &os, const DiGraph &g)
{
    return os << "DiGraph(\n\tforward=" << g.forward.neighbors << ",\n\tbackward=" << g.backward.neighbors << ",\n\tdeleted=" << g.deleted << "\n)";
}

ostream& operator<<(ostream &os, const TwoHopCover &c)
{
    return os << "2HOP(\n\tin=" << c.in << ",\n\tout=" << c.out << "\n)";
}
