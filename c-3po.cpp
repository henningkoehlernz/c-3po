// Implementation of [min,max]-3PO algorithm for constructing 2-hop cover
// compile with g++ -o c3po -O3 -Wall c-3po.cpp

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include <sstream>
#include <string>

#define DEBUG(X) //cout << X << endl

#define MIN_MAX_CENTRALITY
#define UPDATE_CENTRALITY
#define RANDOMIZE

const size_t query_count = 1000000ul;

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

using namespace std;

typedef uint32_t NodeID;

static const uint32_t DEFAULT = UINT32_MAX;
static const uint32_t DELETED = UINT32_MAX - 1;

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
    void prune(const vector<uint32_t> &status)
    {
        if ( lazy_size < vector::size() )
        {
            auto it = begin();
            while ( it != end() )
            {
                if ( status[*it] == DELETED )
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

ostream& operator<<(ostream &os, const LazyList &l);

// one-directional graph
class PartialGraph
{
public:
    vector<LazyList> neighbors;
    vector<NodeID> &last_visited;

    PartialGraph(size_t size, vector<NodeID> &last_visited) : neighbors(size), last_visited(last_visited) {}
};

class DiGraph
{
public:
    PartialGraph forward;
    PartialGraph backward;
    // for tracking deleted & visited nodes during propagation
    vector<NodeID> last_visited;

    DiGraph(size_t size) : forward(size, last_visited), backward(size, last_visited), last_visited(size, DEFAULT) {}

    void insert_edge(NodeID from, NodeID to)
    {
        forward.neighbors[from].push_back(to);
        backward.neighbors[to].push_back(from);
    }

    void remove_node(NodeID node)
    {
        for ( NodeID neighbor : forward.neighbors[node] )
            --backward.neighbors[neighbor];
        for ( NodeID neighbor : backward.neighbors[node] )
            --forward.neighbors[neighbor];
        forward.neighbors[node].clear();
        backward.neighbors[node].clear();
        last_visited[node] = DELETED;
    }

    uint64_t centrality(NodeID node)
    {
        uint64_t in = backward.neighbors[node].size();
        uint64_t out = forward.neighbors[node].size();
#ifdef MIN_MAX_CENTRALITY
        return in <= out ? (in << 32) | out : (out << 32) | in;
#else
        return (in + 1) * (out + 1);
#endif
    }

    size_t size()
    {
        return last_visited.size();
    }
};

ostream& operator<<(ostream &os, const DiGraph &g);

// break centrality-ties using pseudo-random permutation of node ID
class WeightedNode
{
public:
    uint64_t weight;
    uint32_t node;

    WeightedNode(uint64_t weight, uint32_t node) : weight(weight), node(node) {}

    bool operator<(WeightedNode other) const
    {
#ifdef RANDOMIZE
        // "randomized" by priority_queue implementation
        return weight < other.weight;
#else
        return weight < other.weight || (weight == other.weight && node < other.node);
#endif
    }
};

// functions for dealing with label sets
template<typename C>
bool sorted_intersect(const C &ca, const C &cb)
{
    auto a = ca.cbegin();
    auto b = cb.cbegin();
    while ( a != ca.cend() && b != cb.cend() )
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

template<typename C,typename E>
void sorted_insert(C &c, E e)
{
    auto pos = std::lower_bound(c.cbegin(), c.cend(), e);
    c.insert(pos, e);
}

template<typename C,typename E>
void sorted_erase(C &c, E e)
{
    auto pos = std::lower_bound(c.cbegin(), c.cend(), e);
    if ( pos != c.cend() && *pos == e )
        c.erase(pos);
}

typedef vector<NodeID> LabelSet;

bool contains(const LabelSet &labels, NodeID node)
{
    return std::binary_search(labels.cbegin(), labels.cend(), node);
}

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
            if ( contains(in[node], node) )
                ++sum;
            if ( contains(out[node], node) )
                ++sum;
        }
        return sum;
    }

    void map(vector<NodeID> f)
    {
        for ( LabelSet &ls : in )
        {
            for ( NodeID &label : ls )
                label = f[label];
            sort(ls.begin(), ls.end());
        }
        for ( LabelSet &ls : out )
        {
            for ( NodeID &label : ls )
                label = f[label];
            sort(ls.begin(), ls.end());
        }
    }

    // functions for reachability testing

    void erase_self()
    {
        for ( NodeID node = 0; node < in.size(); ++node )
        {
            sorted_erase(in[node], node);
            sorted_erase(out[node], node);
        }
    }

    void insert_self()
    {
        for ( NodeID node = 0; node < in.size(); ++node )
        {
            sorted_insert(in[node], node);
            sorted_insert(out[node], node);
        }
    }

    bool can_reach_all_self(NodeID from, NodeID to) const
    {
        return sorted_intersect(out[from], in[to]);
    }

    bool can_reach_remote(NodeID from, NodeID to) const
    {
        return from == to
            || sorted_intersect(out[from], in[to]);
    }

    bool can_reach_no_self(NodeID from, NodeID to) const
    {
        return from == to
            || std::binary_search(out[from].cbegin(), out[from].cend(), to)
            || std::binary_search(in[to].cbegin(), in[to].cend(), from)
            || sorted_intersect(out[from], in[to]);
    }
};

ostream& operator<<(ostream &os, const TwoHopCover &c);

void propagate_prune(PartialGraph &g, NodeID node, NodeID label, LabelSet &node_labels, vector<LabelSet> &labels)
{
    // propagate remote labels
    vector<NodeID> stack;
    for ( NodeID neighbor : g.neighbors[node] )
        // will be removed after, so just check instead of pruning
        if ( g.last_visited[neighbor] != DELETED )
            stack.push_back(neighbor);
    while ( !stack.empty() )
    {
        NodeID next = stack.back(); stack.pop_back();
        if ( g.last_visited[next] != node && !sorted_intersect(labels[next], node_labels) )
        {
            g.last_visited[next] = node;
            labels[next].push_back(label);
            // prune while we're here
            g.neighbors[next].prune(g.last_visited);
            for ( NodeID neighbor : g.neighbors[next] )
                stack.push_back(neighbor);
        }
    }
    // insert self label
    if ( g.neighbors[node].size() )
        node_labels.push_back(label);
}

TwoHopCover pick_propagate_prune(DiGraph &g, vector<NodeID> &pick_order)
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
        NodeID node = weighted.node;
#ifndef UPDATE_CENTRALITY
        const NodeID label = pick_order.size();
        propagate_prune(g.forward, node, label, labels.out[node], labels.in);
        propagate_prune(g.backward, node, label, labels.in[node], labels.out);
        g.remove_node(node);
        pick_order.push_back(node);
#else
        uint64_t new_weight = g.centrality(node);
        if ( new_weight < weighted.weight )
        {
            weighted.weight = new_weight;
            q.push(weighted);
        }
        else
        {
            DEBUG("picked " << node << " with C=[" << (new_weight >> 32) << ", " << (new_weight & 0xffffffffu) << "]");
            // use ascending labels so that label sets can simply be appended
            const NodeID label = pick_order.size();
            propagate_prune(g.forward, node, label, labels.out[node], labels.in);
            propagate_prune(g.backward, node, label, labels.in[node], labels.out);
            g.remove_node(node);
            pick_order.push_back(node);
        }
#endif
    }
    return labels;
}

DiGraph read_graph(std::istream &in)
{
    size_t nodes;
    string line;
    NodeID node;
    // first line contains #nodes (optionally followed by comment)
    in >> nodes;
    DiGraph g(nodes);
    getline(in, line);
    // remaining lines have the format node child_1 ... child_n
    while ( in ) {
        getline(in, line);
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

DiGraph read_graph_from_file(char* filename)
{
    ifstream ifs(filename);
    return read_graph(ifs);
}

struct Query
{
    NodeID from, to;
    Query(NodeID from, NodeID to) : from(from), to(to) {}
};

Query random_query(size_t nodes)
{
    return Query(std::rand() % nodes, std::rand() % nodes);
}

// convert #labels into MB
double l2mb(size_t labels)
{
    return labels * 4 / (1024.0 * 1024.0);
}

int main (int argc, char *argv[])
{
    DiGraph g = argc > 1 ? read_graph_from_file(argv[1]) : read_graph(std::cin);
    cout << "parsed graph (" << g.size() << " nodes)" << endl;
    if ( g.size() )
    {
        DEBUG("g=" << g);
        vector<NodeID> pick_order;
        pick_order.reserve(g.size());
        auto t_start = chrono::high_resolution_clock::now();
        TwoHopCover cover = pick_propagate_prune(g, pick_order);
        auto t_stop = chrono::high_resolution_clock::now();
        long dur_ms = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        // report stats
        cover.map(pick_order); // needed for identifying self-labels
        DEBUG("cover=" << cover);
        size_t self_labels = cover.self_labels();
        size_t remote_labels = cover.labels() - self_labels;
        cout << "found 2-hop cover with " << remote_labels << " remote and " << self_labels << " self labels in " << static_cast<double>(dur_ms) / 1000.0 << "s" << endl;
        cout << "=> Index size = " << l2mb(remote_labels) << " + " << l2mb(self_labels) << " = " << l2mb(remote_labels + self_labels) << " MB" << endl;
        // test query speed
        vector<Query> queries;
        for ( size_t i = 0; i < query_count; ++i )
            queries.push_back(random_query(g.size()));
        // time all three label storage variants
        t_start = chrono::high_resolution_clock::now();
        for ( Query q : queries )
            cover.can_reach_remote(q.from, q.to);
        t_stop = chrono::high_resolution_clock::now();
        long dur_remote = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "Processed " << query_count << " queries in " << dur_remote << "ms using remote-referenced self-labels" << endl;
        // no self-labels
        cover.erase_self();
        t_start = chrono::high_resolution_clock::now();
        for ( Query q : queries )
            cover.can_reach_no_self(q.from, q.to);
        t_stop = chrono::high_resolution_clock::now();
        long dur_no_self = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "Processed " << query_count << " queries in " << dur_no_self << "ms using no self-labels" << endl;
        // all self-labels
        cover.insert_self();
        t_start = chrono::high_resolution_clock::now();
        for ( Query q : queries )
            cover.can_reach_all_self(q.from, q.to);
        t_stop = chrono::high_resolution_clock::now();
        long dur_all_self = chrono::duration_cast<chrono::milliseconds>(t_stop - t_start).count();
        cout << "Processed " << query_count << " queries in " << dur_all_self << "ms using all self-labels" << endl;
    }
    return 0;
}

// -------------------------- output functions --------------------------------

ostream& operator<<(ostream &os, const LazyList &l)
{
    return os << l.size() << " of " << static_cast<vector<NodeID>>(l);
}

ostream& operator<<(ostream &os, const DiGraph &g)
{
    return os << "DiGraph(\n\tforward=" << g.forward.neighbors << ",\n\tbackward=" << g.backward.neighbors << ",\n\tlast_visited=" << g.last_visited << "\n)";
}

ostream& operator<<(ostream &os, const TwoHopCover &c)
{
    return os << "2HOP(\n\tin=" << c.in << ",\n\tout=" << c.out << "\n)";
}
