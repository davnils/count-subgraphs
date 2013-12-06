#include <boost/graph/graphviz.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <iostream>

#include "TestUtils.hpp"
#include "TreeDecomposition.hpp"

namespace count { namespace test {

/**
 */
static bool isBinaryTree(const tree_decomp_t & tree)
{
  if(boost::num_vertices(tree) == 1)
  {
    return true;
  }

  for(auto it = boost::vertices(tree); it.first != it.second; it.first++)
  {
    auto degree = boost::out_degree(*it.first, tree);
    if(degree < 1 || degree > 3)
    {
      return false;
    }
  }

  return true;
}

/**
 */
void testBinaryDecomposition(std::ostream & os)
{
  const unsigned int TESTS = 1000, SEED = 0xf00d;

  boost::random::mt19937 gen(SEED);

  for(unsigned int test = 0; test < TESTS; ++test)
  {
    os
      << "----------------NEW TEST-------------------\n"
      << "Running test #" << test << "\n";

    auto randomGraph = ::count::test::generateConnectedGraph(gen);
    boost::write_graphviz(os, randomGraph);
    os << std::endl;

    auto decomposed = buildTreeDecomposition(randomGraph);
    ::count::visualizeDecomposition(os, decomposed);
    auto binaryTree = convertToBinaryTree(decomposed);
    //::count::visualizeDecomposition(os, binaryTree.first);

    assert(isBinaryTree(binaryTree.first));

    os
      << "-------------------------------------------\n"
      << std::endl;
  }
}

} }
