#include <fstream>
#include <queue>
#include <sstream>
#include <string>

#include "c-3po.h"

using namespace std;

//------------------------- Neighbor --------------------------------

#ifdef ESTIMATE_ANC_DESC
bool Neighbor::operator<(Neighbor other) {
    return estimate < other.estimate;
}
#endif

Neighbor::Neighbor(NodeID node) : node(node) {
#ifdef ESTIMATE_ANC_DESC
    estimate = DEFAULT;
#endif
}

//------------------------- LazyList --------------------------------

size_t LazyList::size() const
{
    return lazy_size;
}

void LazyList::push_back(NodeID val)
{
    vector::push_back(val);
    ++lazy_size;
}

void LazyList::operator--()
{
    --lazy_size;
}

void LazyList::shrink(const vector<uint32_t> &status, size_t start)
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

void LazyList::clear() // override (non-virtual)
{
    vector::clear();
    lazy_size = 0;
}

// erase non-deleted nodes given as sorted vector
size_t LazyList::erase_all(const vector<NodeID> &sorted)
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

//------------------------- PartialGraph ----------------------------

PartialGraph::PartialGraph(size_t size, vector<NodeID> *last_visited) : neighbors(size), last_visited(last_visited) {}

PartialGraph::PartialGraph(const PartialGraph &g) : neighbors(g.neighbors), last_visited(nullptr) {}

bool PartialGraph::deleted(NodeID node) const {
    return (*last_visited)[node] == DELETED;
}

//------------------------- Estimate --------------------------------

Estimate::Estimate() : estimate(0), tree_estimate(0), tree_sum(0), tree_parent(DEFAULT) {}

//------------------------- DiGraph ---------------------------------

DiGraph::DiGraph(size_t size) : last_visited(size, DEFAULT), forward(size, &last_visited), backward(size, &last_visited) {
#ifdef ESTIMATE_ANC_DESC
    anc_estimate.resize(size);
    desc_estimate.resize(size);
#endif
}

DiGraph::DiGraph(const DiGraph &g) : last_visited(g.last_visited), forward(g.forward), backward(g.backward)
{
#ifdef ESTIMATE_ANC_DESC
    anc_estimate = g.anc_estimate;
    desc_estimate = g.desc_estimate;
#endif
    // fix last_visited pointers
    forward.last_visited = &last_visited;
    backward.last_visited = &last_visited;
}

void DiGraph::insert_edge(NodeID from, NodeID to)
{
    forward.neighbors[from].push_back(to);
    backward.neighbors[to].push_back(from);
}

void DiGraph::remove_node(NodeID node)
{
    for ( Neighbor neighbor : forward.neighbors[node] )
        --backward.neighbors[neighbor.node];
    for ( Neighbor neighbor : backward.neighbors[node] )
        --forward.neighbors[neighbor.node];
    last_visited[node] = DELETED;
}

centrality_t DiGraph::centrality(NodeID node) const
{
#ifdef ESTIMATE_ANC_DESC
    uint64_t in = anc_estimate[node].estimate;
    uint64_t out = desc_estimate[node].estimate;
#else
    uint64_t in = backward.neighbors[node].size();
    uint64_t out = forward.neighbors[node].size();
#endif
#ifdef GAIN_COST_CENTRALITY
    centrality_t gain = in * out + in + out;
    centrality_t cost = in + out + 1e-6;
    return gain / cost;
#else
    return (in + 1) * (out + 1);
#endif
}

size_t DiGraph::size() const
{
    return last_visited.size();
}

#ifdef ESTIMATE_ANC_DESC
void DiGraph::init_estimate_trees()
{
    // requires nodes to be numbered in topological order
    for ( NodeID node = 0; node < size(); ++node )
    {
        LazyList &parents = backward.neighbors[node];
        if ( parents.size() > 0 )
        {
            for ( Neighbor &parent : parents )
                parent.estimate = anc_estimate[parent.node].estimate;
            sort(parents.begin(), parents.end());
            desc_estimate[node].tree_parent = parents.back().node;
            anc_estimate[node].estimate = parents.size() + parents.back().estimate;
        }
        // compute descendants estimates in reverse order
        NodeID rev = size() - node - 1;
        LazyList &children = forward.neighbors[rev];
        if ( children.size() > 0 )
        {
            for ( Neighbor &child : children )
                child.estimate = desc_estimate[child.node].estimate;
            sort(children.begin(), children.end());
            anc_estimate[rev].tree_parent = children.back().node;
            desc_estimate[rev].estimate = children.size() + children.back().estimate;
        }
    }
}

void DiGraph::estimate_anc_desc()
{
    // requires nodes to be numbered in topological order
    for ( NodeID node = 0; node < size(); ++node )
    {
        // init tree estimate
        if ( anc_estimate[node].tree_parent != DEFAULT )
            anc_estimate[anc_estimate[node].tree_parent].tree_sum += anc_estimate[node].tree_sum + 1;
        // compure total estimate
        LazyList &parents = backward.neighbors[node];
        if ( parents.size() > 0 )
        {
            for ( Neighbor &parent : parents )
            {
                parent.estimate = anc_estimate[parent.node].estimate;
                anc_estimate[node].tree_estimate += anc_estimate[parent.node].tree_sum + 1;
            }
            sort(parents.begin(), parents.end());
            anc_estimate[node].estimate = max(static_cast<uint32_t>(parents.size()) + parents.back().estimate, anc_estimate[node].tree_estimate);
        }
        // compute descendants estimates in reverse order
        NodeID rev = size() - node - 1;
        if ( desc_estimate[rev].tree_parent != DEFAULT )
            desc_estimate[desc_estimate[rev].tree_parent].tree_sum += desc_estimate[rev].tree_sum + 1;
        LazyList &children = forward.neighbors[rev];
        if ( children.size() > 0 )
        {
            for ( Neighbor &child : children )
            {
                child.estimate = desc_estimate[child.node].estimate;
                desc_estimate[rev].tree_estimate += desc_estimate[child.node].tree_sum + 1;
            }
            sort(children.begin(), children.end());
            desc_estimate[rev].estimate = max(static_cast<uint32_t>(children.size()) + children.back().estimate, desc_estimate[rev].tree_estimate);
        }
    }
    DEBUG("anc_estimate=" << anc_estimate);
    DEBUG("desc_estimate=" << desc_estimate);
}
#endif

//---------------------- LabelSet utils -----------------------------

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

bool contains(const LabelSet &labels, NodeID node)
{
    return std::binary_search(labels.cbegin(), labels.cend(), node);
}

//------------------------- TwoHopCover -----------------------------

TwoHopCover::TwoHopCover(size_t nodes) : in(nodes), out(nodes) {}

size_t TwoHopCover::labels() const
{
    size_t sum = 0;
    for ( const LabelSet &ls : in )
        sum += ls.size();
    for ( const LabelSet &ls : out )
        sum += ls.size();
    return sum;
}

size_t TwoHopCover::self_labels() const
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

void TwoHopCover::map(vector<NodeID> f)
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

void TwoHopCover::erase_self()
{
    for ( NodeID node = 0; node < in.size(); ++node )
    {
        sorted_erase(in[node], node);
        sorted_erase(out[node], node);
    }
}

void TwoHopCover::insert_self()
{
    for ( NodeID node = 0; node < in.size(); ++node )
    {
        sorted_insert(in[node], node);
        sorted_insert(out[node], node);
    }
}

bool TwoHopCover::can_reach_all_self(NodeID from, NodeID to) const
{
#ifdef TOP_FILTER
    if ( from > to )
        return false;
#endif
    return sorted_intersect(out[from], in[to]);
}

bool TwoHopCover::can_reach_remote(NodeID from, NodeID to) const
{
#ifdef TOP_FILTER
    if ( from > to )
        return false;
#endif
    return from == to
        || sorted_intersect(out[from], in[to]);
}

bool TwoHopCover::can_reach_no_self(NodeID from, NodeID to) const
{
#ifdef TOP_FILTER
    if ( from > to )
        return false;
#endif
    return from == to
        || std::binary_search(out[from].cbegin(), out[from].cend(), to)
        || std::binary_search(in[to].cbegin(), in[to].cend(), from)
        || sorted_intersect(out[from], in[to]);
}

//------------------------- algorithms ------------------------------

void propagate_prune(PartialGraph &g, NodeID node, NodeID label, LabelSet &node_labels, vector<LabelSet> &labels)
{
    // propagate remote labels
    vector<NodeID> stack;
    bool inserted = false;
    g.neighbors[node].shrink(*g.last_visited);
    for ( Neighbor neighbor : g.neighbors[node] )
        stack.push_back(neighbor.node);
    while ( !stack.empty() )
    {
        NodeID next = stack.back(); stack.pop_back();
        if ( (*g.last_visited)[next] != node && !sorted_intersect(labels[next], node_labels) )
        {
            (*g.last_visited)[next] = node;
            labels[next].push_back(label);
            inserted = true;
            g.neighbors[next].shrink(*g.last_visited);
            for ( Neighbor neighbor : g.neighbors[next] )
                stack.push_back(neighbor.node);
        }
    }
    // insert self label
    if ( inserted )
        node_labels.push_back(label);
}

#ifdef ESTIMATE_ANC_DESC
template<class TopCompare>
void update_estimates(PartialGraph &g, PartialGraph &rg, NodeID node, vector<Estimate> &estimates)
{
    priority_queue<NodeID, vector<NodeID>, TopCompare> q;
    for ( Neighbor neighbor : g.neighbors[node] )
        q.push(neighbor.node);
    // update tree estimates (push)
    NodeID tree_anc = node;
    while ( tree_anc != DEFAULT )
    {
        const NodeID tree_parent = estimates[tree_anc].tree_parent;
        if ( tree_parent != DEFAULT )
            estimates[tree_parent].tree_sum -= estimates[node].tree_sum + 1;
        for ( Neighbor neighbor : g.neighbors[tree_anc] )
        {
            // only update estimates of ancestors that might be affected by tree estimate reduction
            if ( estimates[neighbor.node].tree_estimate == estimates[neighbor.node].estimate )
                q.push(neighbor.node);
            estimates[neighbor.node].tree_estimate -= estimates[node].tree_sum + 1;
        }
        tree_anc = tree_parent;
    }
    // make sure node is no longer a tree parent
    for ( Neighbor neighbor : rg.neighbors[node] )
        if ( estimates[neighbor.node].tree_parent == node )
            estimates[neighbor.node].tree_parent = DEFAULT;
    // update final estimates (pull)
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
                uint32_t ne = estimates[rn[pos].node].estimate;
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
                rn.shrink(*g.last_visited, shrink_from);
            sort(rn.begin() + pos, rn.end());
            // ready to re-estimate
            new_estimate = max(static_cast<uint32_t>(rn.size()) + rn.back().estimate, estimates[next].tree_estimate);
        }
        // check if estimated value changed, as transitive edges and tree estimate updates can cause multiple visits
        if ( new_estimate < estimates[next].estimate )
        {
            estimates[next].estimate = new_estimate;
            // shrink again - may not have been visited during label propagation due to label pruning
            g.neighbors[next].shrink(*g.last_visited);
            for ( Neighbor neighbor : g.neighbors[next] )
                q.push(neighbor.node);
        }
    }
}
#endif

// break centrality-ties using pseudo-random permutation of node ID
class WeightedNode
{
public:
    centrality_t weight;
    NodeID node;

    WeightedNode(centrality_t weight, NodeID node) : weight(weight), node(node) {}

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

TwoHopCover pick_propagate_prune(DiGraph &g, vector<NodeID> &pick_order)
{
    TwoHopCover labels(g.size());
#ifdef ESTIMATE_ANC_DESC
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
        centrality_t new_weight = g.centrality(node);
        if ( new_weight < weighted.weight )
        {
            weighted.weight = new_weight;
            q.push(weighted);
        }
        else
        {
            DEBUG("picked " << node << " with C=" << new_weight);
            // use ascending labels so that label sets can simply be appended
            const NodeID label = pick_order.size();
            propagate_prune(g.forward, node, label, labels.out[node], labels.in);
            propagate_prune(g.backward, node, label, labels.in[node], labels.out);
            g.remove_node(node);
#ifdef ESTIMATE_ANC_DESC
            update_estimates<std::greater<NodeID>>(g.forward, g.backward, node, g.anc_estimate);
            update_estimates<std::less<NodeID>>(g.backward, g.forward, node, g.desc_estimate);
            DEBUG("\tanc_estimate=" << g.anc_estimate);
            DEBUG("\tdesc_estimate=" << g.desc_estimate);
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

ostream& operator<<(ostream &os, const Estimate &e)
{
    os << "E(e=" << e.estimate << ", te=" << e.tree_estimate << ", ts=" << e.tree_sum;
    if ( e.tree_parent != DEFAULT )
        os << ", tp=" << e.tree_parent;
    return os << ")";
}

ostream& operator<<(ostream &os, const DiGraph &g)
{
    return os << "DiGraph("
        << "\n\tforward=" << g.forward.neighbors
        << ",\n\tbackward=" << g.backward.neighbors
        << "\n)";
}

ostream& operator<<(ostream &os, const TwoHopCover &c)
{
    return os << "2HOP(\n\tin=" << c.in << ",\n\tout=" << c.out << "\n)";
}
