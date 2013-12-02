#pragma once

#include <boost/graph/adjacency_list.hpp>
#include <set>

namespace count {

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
                              std::set<unsigned int>>
        tree_decomp_t;

typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS>
        undirected_graph_t;

tree_decomp_t buildTreeDecomposition(const undirected_graph_t & inputGraph);
void visualizeDecomposition(std::ostream & os, const tree_decomp_t & tree);

}
