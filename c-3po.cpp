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

// toggles for algorithmic features
#define MIN_MAX_CENTRALITY
#define UPDATE_CENTRALITY
#define RANDOMIZE
#define ESTIMATE_ANC_DESC

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

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2> &p)
{
    return os << "(" << p.first << ", " << p.second << ")";
}

using namespace std;

typedef uint32_t NodeID;

static const uint32_t DEFAULT = UINT32_MAX;
static const uint32_t DELETED = UINT32_MAX - 1;

struct Neighbor
{
    NodeID node;
#ifdef ESTIMATE_ANC_DESC
    uint32_t estimate;
    bool operator<(Neighbor other) { return estimate < other.estimate; }
#endif
    Neighbor(NodeID node) : node(node) {
#ifdef ESTIMATE_ANC_DESC
        estimate = DEFAULT;
#endif
    }
};

// vector of nodes with lazy deletion
class LazyList: public vector<Neighbor>
{
    size_type lazy_size = 0; // track actual data size
public:

    size_type size() const // override (non-virtual)
    {
        return lazy_size;
    }

    void push_back(NodeID val) // override (non-virtual)
    {
        vector::push_back(val);
        ++lazy_size;
    }

    // notify of size reduction
    void operator--()
    {
        --lazy_size;
    }

    // remove deleted nodes while preserving node order (call before iterating over nodes => no complexity increase)
    void shrink(const vector<uint32_t> &status, size_t start = 0)
    {
        if ( lazy_size < vector::size() )
        {
            auto it = begin() + start;
            // find first deleted node - must exist
            while ( status[it->node] != DELETED )
                ++it;
            size_t deleted = 1;
            // shift non-deleted nodes into next free space
            auto next_free = it++;
            while ( it != end() )
            {
                if ( status[it->node] == DELETED )
                {
                    ++it;
                    ++deleted;
                }
                else
                    *(next_free++) = *(it++);
            }
            resize(vector::size() - deleted, 0);
        }
    }

    void clear() // override (non-virtual)
    {
        vector::clear();
        lazy_size = 0;
    }

    // erase non-deleted nodes given as sorted vector
    size_t erase_all(const vector<NodeID> &sorted)
    {
        size_t erased = 0;
        auto it = begin();
        while ( it != end() )
        {
            if ( binary_search(sorted.cbegin(), sorted.cend(), it->node) )
            {
                *it = back();
                pop_back();
                ++erased;
            }
            else
                ++it;
        }
        lazy_size -= erased;
        return erased;
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
    bool deleted(NodeID node) const { return last_visited[node] == DELETED; }
};

class DiGraph
{
public:
    PartialGraph forward;
    PartialGraph backward;
    // for tracking deleted & visited nodes during propagation
    vector<NodeID> last_visited;
#ifdef ESTIMATE_ANC_DESC
    vector<uint32_t> anc_estimate; // total #ancestors
    vector<uint32_t> desc_estimate; // total #descendants
#endif

    DiGraph(size_t size) : forward(size, last_visited), backward(size, last_visited), last_visited(size, DEFAULT) {
#ifdef ESTIMATE_ANC_DESC
        anc_estimate.resize(size, 0);
        desc_estimate.resize(size, 0);
#endif
    }

    void insert_edge(NodeID from, NodeID to)
    {
        forward.neighbors[from].push_back(to);
        backward.neighbors[to].push_back(from);
    }

    void remove_node(NodeID node)
    {
        for ( Neighbor neighbor : forward.neighbors[node] )
            --backward.neighbors[neighbor.node];
        for ( Neighbor neighbor : backward.neighbors[node] )
            --forward.neighbors[neighbor.node];
        last_visited[node] = DELETED;
    }

    uint64_t centrality(NodeID node) const
    {
#ifdef ESTIMATE_ANC_DESC
        uint64_t in = anc_estimate[node];
        uint64_t out = desc_estimate[node];
#else
        uint64_t in = backward.neighbors[node].size();
        uint64_t out = forward.neighbors[node].size();
#endif
#ifdef MIN_MAX_CENTRALITY
        return in <= out ? (in << 32) | out : (out << 32) | in;
#else
        return (in + 1) * (out + 1);
#endif
    }

    size_t size() const
    {
        return last_visited.size();
    }

#ifdef ESTIMATE_ANC_DESC
    void estimate_anc_desc()
    {
        // requires nodes to be numbered in topological order
        for ( NodeID node = 0; node < size(); ++node )
        {
            LazyList &parents = backward.neighbors[node];
            if ( parents.size() > 0 )
            {
                for ( Neighbor &parent : parents )
                    parent.estimate = anc_estimate[parent.node];
                sort(parents.begin(), parents.end());
                anc_estimate[node] = parents.size() + parents.back().estimate;
            }
            // compute descendants estimates in reverse order
            NodeID rev = size() - node - 1;
            LazyList &children = forward.neighbors[rev];
            if ( children.size() > 0 )
            {
                for ( Neighbor &child : children )
                    child.estimate = desc_estimate[child.node];
                sort(children.begin(), children.end());
                desc_estimate[rev] = children.size() + children.back().estimate;
            }
        }
        DEBUG("anc_estimate=" << anc_estimate);
        DEBUG("desc_estimate=" << desc_estimate);
    }
#endif
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
    g.neighbors[node].shrink(g.last_visited);
    for ( Neighbor neighbor : g.neighbors[node] )
        stack.push_back(neighbor.node);
    while ( !stack.empty() )
    {
        NodeID next = stack.back(); stack.pop_back();
        if ( g.last_visited[next] != node && !sorted_intersect(labels[next], node_labels) )
        {
            g.last_visited[next] = node;
            labels[next].push_back(label);
            g.neighbors[next].shrink(g.last_visited);
            for ( Neighbor neighbor : g.neighbors[next] )
                stack.push_back(neighbor.node);
        }
    }
    // insert self label
    if ( g.neighbors[node].size() )
        node_labels.push_back(label);
}

#ifdef ESTIMATE_ANC_DESC
void get_po_tree(const DiGraph &g, vector<vector<NodeID>> &po_tree, vector<NodeID> &parents)
{
    const size_t V = g.size();
    vector<uint32_t> desc_estimate(V);
    po_tree.clear(); po_tree.resize(V + 1); // global root in last element
    parents.clear(); parents.resize(V + 1, V);
    // compute descendants estimates in reverse order
    NodeID node = g.size();
    while ( node-- > 0 )
    {
        const LazyList &children = g.forward.neighbors[node];
        if ( children.size() > 0 )
        {
            NodeID max_node = children[0].node;
            for ( auto child = children.begin() + 1; child != children.end(); ++child )
                if ( desc_estimate[child->node] > desc_estimate[max_node] )
                    max_node = child->node;
            desc_estimate[node] = children.size() + desc_estimate[max_node];
            po_tree[max_node].push_back(node);
            parents[node] = max_node;
        }
        else
            po_tree[g.size()].push_back(node);
    }
}

vector<uint32_t> get_dt_order(const vector<vector<NodeID>> &tree)
{
    const size_t V = tree.size() - 1;
    vector<uint32_t> dt_order(V + 1, 0);
    vector<NodeID> stack;
    for ( auto it = tree[V].rbegin(); it != tree[V].rend(); ++it )
        stack.push_back(*it);
    uint32_t counter = 1;
    while ( !stack.empty() )
    {
        NodeID node = stack.back(); stack.pop_back();
        dt_order[node] = counter++;
        for ( auto it = tree[node].rbegin(); it != tree[node].rend(); ++it )
            stack.push_back(*it);
    }
    return dt_order;
}

// find transitive edges based on buTR algorithm
// see Algorithm 4 in "Accerating reachability query processing based on DAG reduction"
vector<pair<NodeID,NodeID>> get_transitive(const DiGraph &g)
{
    vector<vector<NodeID>> po_tree;
    vector<NodeID> po_parents;
    get_po_tree(g, po_tree, po_parents);
    vector<uint32_t> dt_order = get_dt_order(po_tree);
    DEBUG("po_tree=" << po_tree);
    DEBUG("dt_order" << dt_order);

    const size_t V = g.size();
    vector<uint32_t> flag(V, UINT32_MAX);
    vector<NodeID> edge(V, DEFAULT);

    vector<NodeID> po_stack = po_tree[V];
    vector<NodeID> dfs_stack;
    vector<pair<NodeID,NodeID>> transitive_edges;

    while ( !po_stack.empty() )
    {
        NodeID root = po_stack.back(); po_stack.pop_back();
        NodeID po_parent = po_parents[root];
        DEBUG("root=" << root << ", parent=" << po_parent << ", flag=" << flag);
        uint32_t root_order = dt_order[root];
        uint32_t parent_order = dt_order[po_parent];
        // collect grand-children of root, pruning children in sub-tree
        dfs_stack.clear();
        for ( Neighbor child : g.forward.neighbors[root] )
        {
            if ( flag[child.node] <= parent_order )
            {
                transitive_edges.push_back(make_pair(root, child.node));
                DEBUG("found transitive edge to max sub-tree: " << root << " -> " << child.node);
            }
            else
            {
                edge[child.node] = root;
                flag[child.node] = root_order;
                for ( Neighbor grand_child : g.forward.neighbors[child.node] )
                    dfs_stack.push_back(grand_child.node);
            }
        }
        // DFS over g starting from grand children of root, pruning nodes in max subtree
        while ( !dfs_stack.empty() )
        {
            NodeID node = dfs_stack.back(); dfs_stack.pop_back();
            // ignore nodes already visited or in subtree
            if ( flag[node] > root_order )
            {
                flag[node] = root_order;
                if ( edge[node] == root )
                {
                    transitive_edges.push_back(make_pair(root, node));
                    DEBUG("found transitive edge: " << root << " -> " << node);
                }
                else
                {
                    for ( Neighbor child : g.forward.neighbors[node] )
                        dfs_stack.push_back(child.node);
                }
            }
        }
        // continue DFS on po-tree
        for ( NodeID po_child : po_tree[root] )
            po_stack.push_back(po_child);
    }
    return transitive_edges;
}

size_t remove_transitive(DiGraph &g)
{
    vector<vector<NodeID>> transitive(g.size());
    for ( pair<NodeID,NodeID> edge : get_transitive(g) )
    {
        transitive[edge.first].push_back(edge.second);
        transitive[edge.second].push_back(edge.first);
    }
    // remove transitive neighbors
    size_t erased = 0;
    for ( NodeID node = 0; node < g.size(); node++ )
    {
        vector<NodeID> &t = transitive[node];
        sort(t.begin(), t.end());
        erased += g.forward.neighbors[node].erase_all(t);
        erased += g.backward.neighbors[node].erase_all(t);
    }
    return erased / 2;
}

// construct 2D-index for estimate tree
vector<pair<uint32_t,uint32_t>> index2D(const PartialGraph &g, const PartialGraph &rg)
{
    const size_t V = g.neighbors.size();
    vector<pair<uint32_t,uint32_t>> index(V);
    vector<NodeID> stack;
    for ( NodeID node = 0; node < V; ++node )
        if ( rg.neighbors[node].size() == 0 )
            stack.push_back(node);
    DEBUG("index2D: stack=" << stack);
    uint32_t c_first = 1, c_second = V;
    while ( !stack.empty() )
    {
        NodeID node = stack.back();
        // are we back-tracking?
        if ( index[node].first )
        {
            index[node].second = c_second--;
            stack.pop_back();
        }
        else
        {
            index[node].first = c_first++;
            for ( Neighbor child : g.neighbors[node] )
            {
                const LazyList &rn = rg.neighbors[child.node];
                // check if child in tree, not just in g
                if ( rn.back().node == node )
                    stack.push_back(child.node);
            }
        }
    }
    return index;
}

// find transitive edges that cause double-counting for anc/desc estimates
vector<pair<NodeID,NodeID>> get_tree_transitive(const PartialGraph &g, const vector<pair<uint32_t,uint32_t>> &index)
{
    struct IndexedNode
    {
        pair<uint32_t,uint32_t> label2D;
        NodeID node;
        IndexedNode(pair<uint32_t,uint32_t> label2D, NodeID node) : label2D(label2D), node(node) {}
        bool operator<(const IndexedNode &other) const { return label2D.first < other.label2D.first; }
    };
    vector<IndexedNode> indexed_neighbors;
    // compare neighbors of each node w.r.t. tree-reachability
    vector<pair<NodeID,NodeID>> edges;
    const size_t V = g.neighbors.size();
    for ( NodeID node = 0; node < V; ++node )
    {
        if ( g.neighbors[node].size() > 1 )
        {
            for ( Neighbor neighbor : g.neighbors[node] )
                indexed_neighbors.push_back(IndexedNode(index[neighbor.node], neighbor.node));
            sort(indexed_neighbors.begin(), indexed_neighbors.end());
            uint32_t min_second = indexed_neighbors[0].label2D.second;
            for ( size_t i = 1; i < indexed_neighbors.size(); ++i )
                // sorted by first label => only need to compare second
                if ( min_second < indexed_neighbors[i].label2D.second )
                    edges.push_back(make_pair(node, indexed_neighbors[i].node));
                else
                    min_second = indexed_neighbors[i].label2D.second;
            indexed_neighbors.clear();
        }
    }
    return edges;
}

size_t remove_tree_transitive(DiGraph &g)
{
    vector<pair<uint32_t,uint32_t>> anc_index = index2D(g.forward, g.backward);
    vector<pair<uint32_t,uint32_t>> desc_index = index2D(g.backward, g.forward);
    DEBUG("anc_index=" << anc_index);
    DEBUG("desc_index=" << desc_index);
    // organize transitive edges by node for effient bulk removal
    vector<vector<NodeID>> transitive(g.size());
    for ( pair<NodeID,NodeID> edge : get_tree_transitive(g.forward, anc_index) )
    {
        transitive[edge.first].push_back(edge.second);
        transitive[edge.second].push_back(edge.first);
        DEBUG("found forward t-edge " << edge.first << " -> " << edge.second);
    }
    for ( pair<NodeID,NodeID> edge : get_tree_transitive(g.backward, desc_index) )
    {
        transitive[edge.first].push_back(edge.second);
        transitive[edge.second].push_back(edge.first);
        DEBUG("found backward t-edge " << edge.first << " -> " << edge.second);
    }
    // remove transitive neighbors
    size_t erased = 0;
    for ( NodeID node = 0; node < g.size(); node++ )
    {
        vector<NodeID> &t = transitive[node];
        sort(t.begin(), t.end());
        erased += g.forward.neighbors[node].erase_all(t);
        erased += g.backward.neighbors[node].erase_all(t);
    }
    return erased / 2;
}

template<class TopCompare>
void update_estimates(PartialGraph &g, PartialGraph &rg, NodeID node, vector<uint32_t> &estimates)
{
    priority_queue<NodeID, vector<NodeID>, TopCompare> q;
    for ( Neighbor neighbor : g.neighbors[node] )
        q.push(neighbor.node);
    while ( !q.empty() )
    {
        NodeID next = q.top(); q.pop();
        LazyList &rn = rg.neighbors[next];
        uint32_t new_estimate = 0;
        if ( rn.size() )
        {
            // update estimates for neighbors and re-sort as much as needed
            vector<Neighbor> updated;
            uint32_t min_ne = UINT32_MAX; // minimal neighbor estimate seen so far
            size_t pos = rn.vector::size();
            size_t shrink_from = DEFAULT;
            while ( pos > 0 )
            {
                if ( g.deleted(rn[--pos].node) )
                {
                    shrink_from = pos;
                    continue;
                }
                // update estimates until all preceeding neighbors must have smaller estimates
                uint32_t ne = estimates[rn[pos].node];
                if ( ne < rn[pos].estimate )
                {
                    rn[pos].estimate = ne;
                    if ( ne < min_ne )
                        min_ne = ne;
                }
                else if ( ne <= min_ne )
                {
                    // all preceeding neighbors must have smaller estimates as well
                    break;
                }
            }
            if ( shrink_from != DEFAULT )
                rn.shrink(g.last_visited, shrink_from);
            sort(rn.begin() + pos, rn.end());
            // ready to re-estimate
            new_estimate = rn.size() + rn.back().estimate;
        }
        // check if estimated value changed, as transitive edges can cause multiple visits
        if ( new_estimate < estimates[next] )
        {
            estimates[next] = new_estimate;
            // shrink again - may not have been visited during label propagation due to label pruning
            g.neighbors[next].shrink(g.last_visited);
            for ( Neighbor neighbor : g.neighbors[next] )
                q.push(neighbor.node);
        }
    }
}
#endif

TwoHopCover pick_propagate_prune(DiGraph &g, vector<NodeID> &pick_order)
{
    TwoHopCover labels(g.size());
#ifdef ESTIMATE_ANC_DESC
    //g.estimate_anc_desc(); remove_tree_transitive(g);
    size_t removed = remove_transitive(g);
    cout << "removed " << removed << " transitive edges" << endl;
    g.estimate_anc_desc();
#endif
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
#ifdef ESTIMATE_ANC_DESC
            update_estimates<std::greater<NodeID>>(g.forward, g.backward, node, g.anc_estimate);
            update_estimates<std::less<NodeID>>(g.backward, g.forward, node, g.desc_estimate);
#endif
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
        {
            if ( node_and_children[0] >= node_and_children[child] )
            {
                cerr << "Error: nodes must be numbered in topological order (" << node_and_children[0] << " -> " << node_and_children[child] << ")" << endl;
                exit(0);
            }
            g.insert_edge(node_and_children[0], node_and_children[child]);
        }
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

ostream& operator<<(ostream &os, const Neighbor &n)
{
#ifdef ESTIMATE_ANC_DESC
    if ( n.estimate == DEFAULT )
        return os << n.node;
    return os << "N(" << n.node << ", " << n.estimate << ")";
#else
    return os << n.node;
#endif
}

ostream& operator<<(ostream &os, const LazyList &l)
{
    if ( l.size() < l.vector::size() )
        os << l.size() << " of ";
    return os << static_cast<vector<Neighbor>>(l);
}

ostream& operator<<(ostream &os, const DiGraph &g)
{
    return os << "DiGraph(\n\tforward=" << g.forward.neighbors << ",\n\tbackward=" << g.backward.neighbors << "\n)";
}

ostream& operator<<(ostream &os, const TwoHopCover &c)
{
    return os << "2HOP(\n\tin=" << c.in << ",\n\tout=" << c.out << "\n)";
}
