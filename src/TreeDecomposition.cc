#include <algorithm>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <iostream>
#include <queue>
#include <vector>

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

static std::list<unsigned int> getChildren(const tree_decomp_t & tree,
                                           const unsigned int vertex,
                                           const std::vector<bool> & visited)
{
  std::list<unsigned int> children;
  for(auto it = boost::adjacent_vertices(vertex, tree); it.first != it.second;
  it.first++)
  {
    auto neighbour = *it.first;
    //Ignore any new vertices created after visited
    if(neighbour < visited.size() && !visited.at(neighbour))
    {
      children.push_front(neighbour);
    }
  }
  
  return children;
}

/**
 *
 */
static tree_decomp_t convertToNiceCriteria1(const tree_decomp_t & inputTree,
                                            const unsigned int root)
{
  tree_decomp_t out(inputTree);
  std::queue<unsigned int> work;
  std::vector<bool> visited(boost::num_vertices(inputTree), false);

  work.push(root);
  visited.at(root) = true;
  while(!work.empty())
  {
    auto p = work.front(); work.pop();

    //Split child edges when there are two children available
    auto children = getChildren(out, p, visited);
    if(children.size() == 2)
    {
      auto splitPath = [&out](unsigned int p, unsigned int child)
      {
        auto middle = boost::add_vertex(out[p], out);
        boost::remove_edge(p, child, out);
        boost::add_edge(p, middle, out);
        boost::add_edge(middle, child, out);
      };

      auto q = children.front(); children.pop_front();
      auto r = children.front(); children.pop_front();

      splitPath(p, q);
      splitPath(p, r);

      work.push(q);
      work.push(r);
      visited.at(q) = true;
      visited.at(r) = true;
    }
  }

  return out;
}

/**
 *
 */
static tree_decomp_t convertToNiceCriteria2(const tree_decomp_t & inputTree,
                                            const unsigned int root)
{
  std::cerr << "[DEBUG] CRITERIA 2" << std::endl;
  auto replaceWithPath = [](unsigned int p, unsigned int q, tree_decomp_t & tree)
  {
    auto size = boost::num_vertices(tree);
    std::vector<unsigned int> pqIntersect(size), pqDiff(size), qpDiff(size);

    //Build X_P `ìntersect` X_Q, X_P \ X_Q and X_Q \ X_P
    {
      decltype(pqIntersect.begin()) last;
      auto tp = tree[p];
      auto tq = tree[q];

      last = std::set_intersection(tp.begin(), tp.end(), tq.begin(), tq.end(),
                                   pqIntersect.begin());
      pqIntersect.resize(last - pqIntersect.begin());

      last = std::set_difference(tp.begin(), tp.end(), tq.begin(), tq.end(),
                                 pqDiff.begin());
      pqDiff.resize(last - pqDiff.begin());

      last = std::set_difference(tq.begin(), tq.end(), tp.begin(), tp.end(),
                                 qpDiff.begin());
      qpDiff.resize(last - qpDiff.begin());
    }

    //Extend p-q edge with path from definition i.e. p - extension - q
    //Build extension from the middle pq vertex and outwards
    auto pq = boost::add_vertex(std::set<unsigned int>(pqIntersect.begin(),
                                pqIntersect.end()), tree);
    auto buildExtension = [&tree, pq]
      (const std::vector<unsigned int> & addition, unsigned int finalVertex)
    {
      auto previous = pq;
      for(auto v : addition)
      {
        auto newBag = tree[previous];
        newBag.insert(v);
        auto ext = boost::add_vertex(newBag, tree);
        boost::add_edge(ext, previous, tree);
        previous = ext;
      }

      boost::add_edge(previous, finalVertex, tree);
    };

    boost::remove_edge(p, q, tree);
    buildExtension(pqDiff, p);
    buildExtension(qpDiff, q);
  };

  auto contractEdge = [](unsigned int p, unsigned int q, tree_decomp_t & tree)
  {
    for(auto it = boost::adjacent_vertices(q, tree); it.first != it.second;
        ++it.first)
    {
      auto next = *it.first;
      if(next != p)
      {
        boost::remove_edge(q, next, tree);
        boost::add_edge(p, next, tree);
      }
    }
  };

  tree_decomp_t out(inputTree);
  std::queue<unsigned int> work;
  std::vector<bool> visited(boost::num_vertices(inputTree), false);

  work.push(root);
  visited.at(root) = true;
  while(!work.empty())
  {
    auto p = work.front(); work.pop();

    //Process all nodes having one children
    auto children = getChildren(out, p, visited);
    if(children.size() != 1)
    {
      //Push all children and continue search
      while(!children.empty())
      {
        work.push(children.front());
        visited.at(children.front()) = true;
        children.pop_front();
      }
      continue;
    }

    //Check if bags are equal, two different cases
    auto q = children.front();
    if(out[p] == out[q])
    {
      //Contract p-q edge and reconsider the updated p vertex
      contractEdge(p, q, out);
      work.push(p);
    }
    else
    {
      replaceWithPath(p, q, out);
      work.push(q);
      visited.at(q) = true;
    }
  }

  //Filter out disconnected vertices
  class VertexFilter
  {
  public:
    VertexFilter() = default;
    VertexFilter(tree_decomp_t tree) : m_tree(tree)
    {
    }

    bool operator()(const unsigned int vertex) const
    {
      return boost::out_degree(vertex, m_tree) > 0;
    }

  private:
    tree_decomp_t m_tree;
  };

  VertexFilter filter(out);
  boost::filtered_graph<tree_decomp_t, boost::keep_all, VertexFilter>
    copy(out, boost::keep_all(), filter);
  tree_decomp_t final;
  boost::copy_graph(copy, final);
  return final;
}

/**
 *
 */
static tree_decomp_t convertToCriteria3(const tree_decomp_t & inputTree,
                                        const unsigned int root)
{
  std::cerr << "[DEBUG] CRITERIA 3" << std::endl;
  tree_decomp_t out(inputTree);
  std::queue<unsigned int> work;
  std::vector<bool> visited(boost::num_vertices(inputTree), false);

  work.push(root);
  visited.at(root) = true;
  while(!work.empty())
  {
    auto p = work.front(); work.pop();
    auto children = getChildren(out, p, visited);

    if(children.empty())
    {
      //Extend child vertex with a path
      //containing decreasing bag sizes
      unsigned int previous = p;
      while(out[previous].size() > 1)
      {
        auto bag = out[previous];
        bag.erase(bag.begin());
        auto next = boost::add_vertex(bag, out);
        boost::add_edge(previous, next, out);
        previous = next;
      }
    }

    while(!children.empty())
    {
      work.push(children.front());
      visited.at(children.front()) = true;
      children.pop_front();
    }
  }

  return out;
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
  auto tree = convertToBinaryTree(inputTree);
  auto niceDecomp = convertToCriteria3(
                      convertToNiceCriteria2(
                        convertToNiceCriteria1(tree.first, tree.second),
                        tree.second),
                      tree.second);
  return std::make_pair(niceDecomp, tree.second);
}

void visualizeDecomposition(std::ostream & os, const tree_decomp_t & tree)
{
  boost::write_graphviz(os, tree, VertexSetWriter(tree));
}

}
