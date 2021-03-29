// compile with g++ -o test -O3 -Wall test.cpp c-3po.cpp tr.cpp

#include <cstdlib>
#include <set>

#include "c-3po.h"
#include "tr.h"

using namespace std;

DiGraph random_dag(size_t nodes, size_t edges)
{
    DiGraph g(nodes);
    // track neighbors as set to avoid duplicates
    vector<set<NodeID>> children(nodes);
    for ( size_t i = 0; i < edges; i++ )
    {
        NodeID from = rand() % nodes;
        NodeID to = rand() % nodes;
        if ( from > to )
            swap(from, to);
        if ( from < to && children[from].find(to) == children[from].end() )
        {
            children[from].insert(to);
            g.insert_edge(from, to);
        }
    }
    return g;
}

vector<DiGraph> random_dags(size_t repeats=10)
{
    vector<size_t> node_nrs { 3, 4, 5, 7, 9, 12, 15, 20 };
    vector<double> densities { 1.0, 1.2, 1.5, 2.0, 3.0 };
    vector<DiGraph> dags;
    for ( size_t node_nr : node_nrs )
        for ( double density : densities )
            for ( size_t i = 0; i < repeats; i++ )
                dags.push_back(random_dag(node_nr, node_nr * density));
    return dags;
}

bool can_reach_dfs(const DiGraph &g, NodeID from, NodeID to)
{
    vector<bool> visited(g.size(), false);
    vector<NodeID> stack { from };
    while ( !stack.empty() )
    {
        NodeID node = stack.back(); stack.pop_back();
        if ( node == to )
            return true;
        if ( !visited[node] )
        {
            visited[node] = true;
            for ( Neighbor neighbor : g.forward.neighbors[node] )
                stack.push_back(neighbor.node);
        }
    }
    return false;
}

vector<pair<NodeID,NodeID>> node_pairs(size_t nodes)
{
    vector<pair<NodeID, NodeID>> pairs;
    for ( NodeID first = 0; first < nodes; ++first )
        for ( NodeID second = 0; second < nodes; ++second )
        pairs.push_back(make_pair(first, second));
    return pairs;
}

void test_ppp(DiGraph &g)
{
    // check actual reachability before graph is consumed
    vector<pair<NodeID,NodeID>> queries = node_pairs(g.size());
    vector<bool> can_reach;
    for ( pair<NodeID,NodeID> query : queries )
        can_reach.push_back(can_reach_dfs(g, query.first, query.second));
    // build 2-hop index
    vector<NodeID> pick_order;
    TwoHopCover cover = pick_propagate_prune(g, pick_order);
    // compare
    for ( size_t q = 0; q < queries.size(); ++q )
    {
        bool cover_reachable = cover.can_reach_remote(queries[q].first, queries[q].second);
        if ( cover_reachable != can_reach[q] )
        {
            cout << "For g=" << g << "\ngot cover=" << cover << "\nwith can_reach_remote(" << queries[q].first << "," << queries[q].second << ")=" << cover_reachable << endl;
            exit(0);
        }
    }
}

int main (int argc, char *argv[])
{
    cout << "Testing pick-propagate-prune ..." << endl;
    for ( DiGraph &g : random_dags() )
        test_ppp(g);
    return 0;
}
