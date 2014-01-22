#pragma once

#include <boost/random.hpp>

#include "TreeDecomposition.hpp"

namespace Count { namespace Test {

Tree::undirected_graph_t generateConnectedGraph(
  boost::random::mt19937 &,
  unsigned int
  );

void runTests(
  std::function<void(void)>,
  const std::string &,
  std::ostream &
  );

} }
