#pragma once

#include "TreeDecomposition.hpp"
#include "Utils.hpp"

namespace Count { namespace Homomorphism {

typedef unsigned long long count_t;
typedef std::map<unsigned int, unsigned int> map_t;

std::list<unsigned int> calculateStingyOrdering(
  const Tree::tree_decomp_t &,
  const unsigned int);

unsigned long long countHomomorphisms(
  const Tree::undirected_graph_t &,
  const Tree::tree_decomp_t &,
  const unsigned int,
  const Tree::undirected_graph_t &,
  const map_t &);

bool isValidPartialHomomorphism(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &,
  const map_t &);

} }
