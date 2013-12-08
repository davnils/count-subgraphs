#pragma once

#include "TreeDecomposition.hpp"

namespace homomorphism {

typedef unsigned long long count_t;

unsigned long long countHomomorphisms(const count::undirected_graph_t &,
                                      const count::tree_decomp_t &,
                                      const count::undirected_graph_t &);

}
