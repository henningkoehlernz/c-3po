// functions for computing transitive reduction
#pragma once

#include "c-3po.h"

void get_po_tree(const DiGraph &g, vector<vector<NodeID>> &po_tree, vector<NodeID> &parents);
vector<uint32_t> get_dt_order(const vector<vector<NodeID>> &tree);

// find transitive edges based on buTR algorithm
// see Algorithm 4 in "Accerating reachability query processing based on DAG reduction"
vector<pair<NodeID,NodeID>> get_transitive(const DiGraph &g);

size_t remove_transitive(DiGraph &g);

// construct 2D-index for estimate tree
vector<pair<uint32_t,uint32_t>> index2D(const PartialGraph &g, const PartialGraph &rg);

// find transitive edges that cause double-counting for anc/desc estimates
vector<pair<NodeID,NodeID>> get_tree_transitive(const PartialGraph &g, const vector<pair<uint32_t,uint32_t>> &index);

size_t remove_tree_transitive(DiGraph &g);
