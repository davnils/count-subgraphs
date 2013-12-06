#include <boost/graph/graphviz.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <queue>

#include "TreeDecomposition.hpp"

namespace count {

class VertexSetWriter
{
public:
  VertexSetWriter(const tree_decomp_t & t)
    : m_tree(t)
  {
  }

  template <typename Item>
  void operator()(std::ostream & os, const Item & v) const
  {
    os << "[label=\"{";
    bool first = true;
    for(auto u : m_tree[v])
    {
      if(!first)
      {
        os << ",";
      }
      first = false;

      os << u;
    }
    os <<  "}\"]";
  }

private:
  const tree_decomp_t & m_tree;
};

tree_decomp_t buildTreeDecomposition(const undirected_graph_t & inputGraph)
{
  dlib::graph<int>::kernel_1a_c inputTransformed;
  inputTransformed.set_number_of_nodes(boost::num_vertices(inputGraph));

  for(auto it = boost::edges(inputGraph); it.first != it.second; it.first++)
  {
    inputTransformed.add_edge(boost::source(*it.first, inputGraph),
                              boost::target(*it.first, inputGraph));
  }

  typedef dlib::set<unsigned long>::compare_1b_c set_t;
  dlib::graph<set_t, set_t>::kernel_1a_c internalDecomp;

  dlib::create_join_tree(inputTransformed, internalDecomp);

  tree_decomp_t outputDecomp(internalDecomp.number_of_nodes());
  for(auto v = 0; v < internalDecomp.number_of_nodes(); v++)
  {
    auto & it = internalDecomp.node(v).data;
    it.reset();
    while(it.move_next())
    {
      outputDecomp[v].insert(it.element());
    }

    for(auto u = v; u < internalDecomp.number_of_nodes(); u++)
    {
      if(internalDecomp.has_edge(v, u))
      {
        boost::add_edge(v, u, outputDecomp);
      }
    }
  }

  return outputDecomp;
}


/**
 * .
 */
static void mergeBranches(tree_decomp_t & tree, unsigned int vertex,
                          unsigned int move, unsigned int target)
{
  std::vector<unsigned int> bagUnion(tree[move].size() + tree[target].size());;
  auto last = std::set_union(tree[move].begin(), tree[move].end(),
                             tree[target].begin(), tree[target].end(),
                             bagUnion.begin());

  tree[target] = std::set<unsigned int>(bagUnion.begin(), last);
  boost::remove_edge(vertex, move, tree);
  boost::add_edge(target, move, tree);
}

/**
 * Convert an existing tree decomposition into a binary tree.
 * Traverse the tree once and push additional subtrees
 * onto neighbouring branches, while maintaining invariants.
 */
std::pair<tree_decomp_t, unsigned int> convertToBinaryTree(const tree_decomp_t & inputTree)
{
  //Handle isolated vertices
  if(boost::num_vertices(inputTree) == 1)
  {
    return std::make_pair(inputTree, 0);
  }

  //Locate some vertex with degree 1
  unsigned int startVertex;
  for(startVertex = 0; startVertex < boost::num_vertices(inputTree); ++startVertex)
  {
    if(boost::out_degree(startVertex, inputTree) == 1)
    {
      break;
    }
  }
  assert(startVertex != boost::num_vertices(inputTree));

  //Traverse the tree and gradually build the corresponding binary tree.
  //Transform the input tree by moving edges and taking union between sets.
  std::vector<bool> visited(boost::num_vertices(inputTree), false);
  tree_decomp_t tree(inputTree);
  std::queue<unsigned int> work;

  work.push(startVertex);
  visited[startVertex] = true;
  while(!work.empty())
  {
    auto vertex = work.front(); work.pop();

    std::list<unsigned int> neighbours;
    for(auto it = boost::adjacent_vertices(vertex, tree); it.first != it.second;
        ++it.first)
    {
      if(!visited.at(*it.first))
      {
        neighbours.push_front(*it.first);
      }
    }

    //Skip vertices satisfying the criteria
    if(boost::out_degree(vertex, tree) <= 3)
    {
      for(auto v : neighbours)
      {
        visited.at(v) = true;
        work.push(v);
      }
      continue;
    }

    //Must be the case that at least three branches,
    //since there is at most once 'visited' vertex.
    assert(neighbours.size() >= 3);

    //Finally move all but two branches
    auto mergeTarget = neighbours.front(); neighbours.pop_front();
    auto secondChild = neighbours.front(); neighbours.pop_front();

    work.push(mergeTarget);
    visited.at(mergeTarget) = true;
    work.push(secondChild);
    visited.at(secondChild) = true;

    for(auto neighbour : neighbours)
    {
      mergeBranches(tree, vertex, neighbour, mergeTarget);
    }
  }

  return std::make_pair(tree, startVertex);
}

/**
 * Converts an existing valid tree decomposition into a nice tree decomposition.
 * A nice tree decomposition satisifies (Bodlaender and Kloks):
 * (1) Every node of the tree has at most two children
 * (2) if a node i has two children j; h then X_i = X_j = X_h,
 * (3) if a node i has one child, then either |X_i|=|X_j| + 1
 *     and X_j \subset X_i or |X_i|=|X_j| − 1 and X_i \subset X_j.
 */
std::pair<tree_decomp_t, unsigned int> convertToNiceDecomposition(const tree_decomp_t & inputTree)
{
  auto binaryTree = convertToBinaryTree(inputTree);
  //TODO
  return binaryTree;
}

void visualizeDecomposition(std::ostream & os, const tree_decomp_t & tree)
{
  boost::write_graphviz(os, tree, VertexSetWriter(tree));
}

}
