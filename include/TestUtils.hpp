#pragma once

#include <boost/random.hpp>

#include "TreeDecomposition.hpp"

namespace Count { namespace Test {

Tree::undirected_graph_t generateConnectedGraph(boost::random::mt19937 &);

} }
