#pragma once

#include <boost/random.hpp>

#include "TreeDecomposition.hpp"

namespace count { namespace test {

undirected_graph_t generateConnectedGraph(boost::random::mt19937 &);

} }
