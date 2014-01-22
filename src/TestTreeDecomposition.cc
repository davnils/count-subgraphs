#include <algorithm>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <functional>
#include <iostream>
#include <list>
#include <queue>
#include <vector>

#include "TestUtils.hpp"
#include "TreeDecomposition.hpp"

using namespace Count::Test;

namespace Count { namespace Tree { namespace Test { 

/**
 * Check if the provided tree is binary.
 *
 * @param tree Input tree.
 * @return True if the tree is valid.
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
 * Check if the provided decomposition is nice.
 *
 * @param tree Input tree.
 * @param root Root vertex.
 * @return True if the tree is valid.
 */
static bool isNiceTreeDecomposition(const tree_decomp_t & tree, unsigned int root)
{
  //check invariant 1
  if(!isBinaryTree(tree))
  {
    return false;
  }

  //check invariants 2 and 3 by performing a BFS traversal from the root
  std::vector<bool> visited(boost::num_vertices(tree), false);
  std::queue<unsigned int> work;

  work.push(root);
  visited.at(root) = true;

  while(!work.empty())
  {
    auto vertex = work.front(); work.pop();
    auto it = boost::adjacent_vertices(vertex, tree);

    std::list<unsigned int> childs(it.first, it.second);
    for(auto & u : childs)
    {
      if(visited.at(u))
      {
        childs.remove(u);
        break;
      }
    }

    //one child case
    if(childs.size() == 1)
    {
      auto const & c = tree[childs.front()];
      auto const & v = tree[vertex];

      auto first  = v.size() == c.size() - 1 &&
                    std::includes(c.begin(), c.end(), v.begin(), v.end());
      auto second = v.size() == c.size() + 1 &&
                    std::includes(v.begin(), v.end(), c.begin(), c.end());

      if(!first && !second)
      {
        return false;
      }
      work.push(childs.front());
      visited.at(childs.front()) = true;
    }
    //two children case
    else if(childs.size() == 2)
    {
      auto first = childs.front(); childs.pop_front();
      auto second = childs.front(); childs.pop_front();
      if(tree[first] != tree[second] || tree[first] != tree[vertex])
      {
        return false;
      }

      work.push(first);
      work.push(second);
      visited.at(first) = true;
      visited.at(second) = true;
    }
  }

  return true;
}

/**
 * Verify that the processed tree decompositions are binary.
 *
 * @param os Output stream to be written.
 */
void testBinaryDecomposition(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen]()
  {
    auto randomGraph = ::Count::Test::generateConnectedGraph(gen, 200);
    auto decomposed = buildTreeDecomposition(randomGraph);
    auto binaryTree = convertToBinaryTree(decomposed);

    assert(isBinaryTree(binaryTree.first));
  };

  runTests(testCase, "BinaryDecomposition", os);
}

/**
 * Verify that the processed tree decompositions are nice.
 *
 * @param os Output stream to be written.
 */
void testNiceTreeDecomposition(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    auto randomGraph = ::Count::Test::generateConnectedGraph(gen, 200);
    auto decomp = buildTreeDecomposition(randomGraph);
    auto niceTree = convertToNiceDecomposition(decomp);

    assert(isBinaryTree(niceTree.first));
    assert(niceTree.first[niceTree.second].empty());
    assert(isNiceTreeDecomposition(niceTree.first, niceTree.second));
  };

  runTests(testCase, "NiceTreeDecomposition", os);
}

} } }
