#pragma once

#include <map>

#include "TreeDecomposition.hpp"
#include "Utils.hpp"

namespace Count { namespace Naive {

unsigned long countHomomorphisms(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &
  );

unsigned long countInjective(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &
  );

unsigned long evaluateDisjointTriple(
  const std::pair<Utils::inj_homo_t, Utils::partition_triple_t> &,
  const Tree::undirected_graph_t &,
  const Utils::PatternPartition &
  );

unsigned long evaluateDisjointTriples(
  const std::map<Utils::inj_homo_t, Utils::partition_triple_t> &,
  const Tree::undirected_graph_t &,
  const Utils::PatternPartition &
  );

} }
