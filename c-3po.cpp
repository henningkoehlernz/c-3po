// Implementation of [min,max]-3PO algorithm for constructing 2-hop cover
// compile with g++ -o c3po -O3 -Wall c-3po.cpp

#include <vector>
#include <cstdint>

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
};

// one-directional graph
class PartialGraph
{
public:
    vector<LazyList> neighbors;
    const vector<bool> &deleted;

    PartialGraph(vector<bool> &deleted) : deleted(deleted) {}
};

class DiGraph
{
public:
    PartialGraph forward;
    PartialGraph backward;
    vector<bool> deleted;

    DiGraph() : forward(deleted), backward(deleted) {}

    void remove_node(NodeID node)
    {
        for ( NodeID neighbor : forward.neighbors[node] )
            --backward.neighbors[neighbor];
        for ( NodeID neighbor : backward.neighbors[node] )
            --forward.neighbors[neighbor];
        deleted[node] = true;
    }

    uint64_t centrality(NodeID node)
    {
        uint64_t in = backward.neighbors[node].size();
        uint64_t out = forward.neighbors[node].size();
        return in <= out ? (in << 32) & out : (out << 32) & in;
    }

};


int main (int argc, char *argv[])
{
    return 0;
}
