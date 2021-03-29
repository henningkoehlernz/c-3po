// functions for computing transitive reduction
#pragma once

#include "c-3po.h"

void get_po_tree(const DiGraph &g, std::vector<std::vector<NodeID>> &po_tree, std::vector<NodeID> &parents);
std::vector<uint32_t> get_dt_order(const std::vector<std::vector<NodeID>> &tree);

// find transitive edges based on buTR algorithm
// see Algorithm 4 in "Accerating reachability query processing based on DAG reduction"
std::vector<std::pair<NodeID,NodeID>> get_transitive(const DiGraph &g);

size_t remove_transitive(DiGraph &g);

// construct 2D-index for estimate tree
std::vector<std::pair<uint32_t,uint32_t>> index2D(const PartialGraph &g, const PartialGraph &rg);

// find transitive edges that cause double-counting for anc/desc estimates
std::vector<std::pair<NodeID,NodeID>> get_tree_transitive(const PartialGraph &g, const std::vector<std::pair<uint32_t,uint32_t>> &index);

size_t remove_tree_transitive(DiGraph &g);
