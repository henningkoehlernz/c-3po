// Run [min,max]-3PO algorithm for constructing 2-hop cover from command line
// compile with g++ -o c3po -O3 -Wall main.cpp c-3po.cpp tr.cpp
// configure flags in c-3po.h

#include <chrono>
#include "c-3po.h"
#include "tr.h"

using namespace std;

const size_t query_count = 1000000ul;

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
#ifdef ESTIMATE_ANC_DESC
        //g.estimate_anc_desc(); remove_tree_transitive(g);
        size_t removed = remove_transitive(g);
        cout << "removed " << removed << " transitive edges" << endl;
#endif
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
