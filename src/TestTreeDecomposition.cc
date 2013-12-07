#include <algorithm>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <iostream>
#include <list>
#include <queue>
#include <vector>

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
static bool isNiceTreeDecomposition(const tree_decomp_t & tree, unsigned int root)
{
  //Check invariant 1
  if(!isBinaryTree(tree))
  {
    return false;
  }

  //Check invariants 2 and 3 by performing a BFS traversal from the root
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

    //One child case
    if(childs.size() == 1)
    {
      auto const & c = tree[childs.front()];
      auto const & v = tree[vertex];

      auto first  = v.size() == c.size() + 1 &&
                    std::includes(c.begin(), c.end(), v.begin(), v.end());
      auto second = v.size() == c.size() - 1 &&
                    std::includes(v.begin(), v.end(), c.begin(), c.end());

      if(!first && !second)
      {
        return false;
      }
      work.push(childs.front());
    }
    //Two children case
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
