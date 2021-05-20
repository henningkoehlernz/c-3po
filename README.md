# Minimal acyclic 2-hop covers for large graphs

Algorithm to construct a 2-hop cover for large sparse directed acyclic graphs, allowing reachability queries to be answered quickly.

## Submodules

To install submodules containting DAGs for testing, use the following commands:

    git submodule init
    git submodule update

Note that this will download over 1 GB worth of graph data.

## Running the code

Run `make` to build, then `./c3po < sample.dag` or `./c3po sample.dag` to run the algorithm with sample.dag as input.

To run the algorithm for all graphs in the submodule, use `sh run.sh`.

To run variants of the algorithm, edit c-3po.h and (un)comment the relevant `#define` statements.
