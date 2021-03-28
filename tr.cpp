#include "tr.h"

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
