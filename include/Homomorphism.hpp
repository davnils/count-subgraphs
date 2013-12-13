#pragma once

#include "TreeDecomposition.hpp"

namespace Homomorphism {

typedef unsigned long long count_t;

unsigned long long countHomomorphisms(const Tree::undirected_graph_t &,
                                      const Tree::tree_decomp_t &,
                                      const unsigned int,
                                      const Tree::undirected_graph_t &);

}
