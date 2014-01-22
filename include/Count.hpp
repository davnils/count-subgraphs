#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <set>

#include "TreeDecomposition.hpp"
#include "Utils.hpp"

namespace Count {

struct subgraph_result
{
  subgraph_result(const Utils::PatternPartition & p) : part(p) {}
  Utils::PatternPartition part;
  std::map<Utils::inj_homo_t, Utils::partition_triple_t> counts;
};

std::map<std::set<unsigned int>, unsigned long> countInjective(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &,
  const std::set<unsigned int> &,
  const Utils::inj_homo_t &
  );


subgraph_result countSubgraphTriples(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &
  );

subgraph_result countSubgraphTriples(
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &,
  const Utils::PatternPartition &
  );

std::vector<int> solveLinearSystem(
  const boost::numeric::ublas::matrix<int> &,
  const boost::numeric::ublas::matrix<int> &
  );

long getCoefficient(
  const unsigned int,
  const unsigned int,
  const unsigned int
  );

std::pair<bool, unsigned long long> countIsoSubgraphs (
  const Tree::undirected_graph_t &,
  const Tree::undirected_graph_t &
  );

}
