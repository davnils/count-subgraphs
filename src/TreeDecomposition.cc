#include <algorithm>
#include <boost/graph/copy.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graphviz.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <iostream>
#include <queue>
#include <vector>

#include "InducedGraph.hpp"
#include "TreeDecomposition.hpp"

namespace Count { namespace Tree {

/**
 * Decompose a component of an undirected graph.
 *
 * @param inputGraph Graph containing a single component.
 * @return The resulting tree decomposition.
 */
static tree_decomp_t decomposeComponent(const undirected_graph_t & inputGraph)
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
  for(auto v = 0u; v < internalDecomp.number_of_nodes(); v++)
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
 * Build a tree decomposition of an undirected simple graph.
 * Multiple components are handled by joining the decompositions.
 *
 * @param inputGraph Graph to be decomposed.
 * @return The resulting tree decomposition.
 */
tree_decomp_t buildTreeDecomposition(const undirected_graph_t & inputGraph)
{
  std::vector<unsigned int> indexToComponent(boost::num_vertices(inputGraph));
  auto components = boost::connected_components(inputGraph, &indexToComponent[0]);

  std::vector<std::set<unsigned int>> componentSubgraphs(components);
  for(unsigned int i = 0; i < boost::num_vertices(inputGraph); ++i)
  {
    componentSubgraphs.at(indexToComponent.at(i)).insert(i);
  }

  //decompose every component separately
  std::vector<tree_decomp_t> decompositions;
  for(auto const & vertices : componentSubgraphs)
  {
    auto subgraph = InducedGraph::buildInducedSubGraphTagged(inputGraph, vertices);

    //build map from component onto original graph vertices
    std::vector<int> vertexMap(boost::num_vertices(subgraph.first), -1);
    for(unsigned int v = 0; v < subgraph.second.size(); ++v)
    {
      if(subgraph.second[v] != -1)
      {
        vertexMap.at(subgraph.second[v]) = v;
      }
    }

    auto decomp = decomposeComponent(subgraph.first);

    //traverse and update every bag to use 'inputGraph' indices
    for(unsigned int v = 0; v < boost::num_vertices(decomp); ++v)
    {
      std::set<unsigned int> newBag;
      for(auto bagItem : decomp[v])
      {
        auto oldVertex = vertexMap.at(bagItem);
        assert(oldVertex != -1);
        newBag.insert(oldVertex);
      }
      decomp[v] = newBag;
    }

    decompositions.push_back(decomp);
  }

  //locates a zero- or one-degree vertex
  auto findDegreeOneVertex =
    []
    (const tree_decomp_t & tree)
  {
    for(auto v = 0u; v < boost::num_vertices(tree); ++v)
    {
      if(boost::out_degree(v, tree) <= 1)
      {
        return v;
      }
    }

    assert(false);
  };

  //merge the decompositions
  while(decompositions.size() != 1)
  {
    auto & outGraph = decompositions.front();
    auto & source = decompositions.back();

    auto oneVertex1 = findDegreeOneVertex(source);
    auto oneVertex2 = findDegreeOneVertex(outGraph);

    //update bags and edges accordingly
    auto offset = boost::num_vertices(outGraph);
    for(auto v = 0u; v < boost::num_vertices(source); ++v)
    {
      boost::add_vertex(source[v], outGraph);
    }

    for(auto it = boost::edges(source); it.first != it.second; ++it.first)
    {
      boost::add_edge(offset + boost::source(*it.first, source),
                      offset + boost::target(*it.first, source),
                      outGraph);
    }

    auto link = boost::add_vertex(std::set<unsigned int>(), outGraph);
    boost::add_edge(link, oneVertex1 + offset, outGraph);
    boost::add_edge(link, oneVertex2, outGraph);

    decompositions.pop_back();
  }

  return decompositions.at(0);
}


/**
 * Moves a subtree based on a source edge to some other vertex.
 *
 * @param tree Tree to be processed.
 * @param vertex Vertex originally connected with the subtree.
 * @param move First vertex in the subtree.
 * @param target New vertex being the parent to the subtree.
 */
static void mergeBranches(
  tree_decomp_t & tree,
  unsigned int vertex,
  unsigned int move,
  unsigned int target
  )
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
 *
 * @param inputTree Tree to be processed.
 * @return The resulting binary decomposition paired with a root vertex.
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
 * Retrieve the children of some vertex in the supplied tree.
 * Will ignore any vertices outside of the visited array, or 
 * any vertices that are marked as visited.
 *
 * @param tree Tree to be queried.
 * @param vertex Parent vertex.
 * @param visited Array indexed by vertices,
 *                indicating if some vertex has been visited.
 */
static std::list<unsigned int> getChildren(
  const tree_decomp_t & tree,
  const unsigned int vertex,
  const std::vector<bool> & visited
  )
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
 * Process the tree to satisfy critieria 1 of nice tree decompositions.
 * 
 * @param inputTree Tree to be processed.
 * @param root Root node of the tree.
 * @return Processed tree satisfying critieria 1.
 */
static tree_decomp_t convertToNiceCriteria1(
  const tree_decomp_t & inputTree,
  const unsigned int root
  )
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
    else if(children.size() == 1)
    {
      work.push(children.front());
      visited.at(children.front()) = true;
    }
  }

  return out;
}

/**
 * Process the tree to satisfy critieria 2 of nice tree decompositions.
 * 
 * @param inputTree Tree to be processed.
 * @param root Root node of the tree.
 * @return Processed tree satisfying critieria 2.
 */
static tree_decomp_t convertToNiceCriteria2(
  const tree_decomp_t & inputTree,
  const unsigned int root
  )
{
  unsigned long largestBag = 0;
  for(auto v = 0u; v < boost::num_vertices(inputTree); ++v)
  {
    largestBag = std::max(largestBag, inputTree[v].size());
  }

  auto replaceWithPath =
    [&largestBag]
    (unsigned int p, unsigned int q, tree_decomp_t & tree)
  {
    auto size = largestBag;
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

    unsigned int pq;
    auto pqBag = std::set<unsigned int>(pqIntersect.begin(), pqIntersect.end());
    if(pqBag == tree[p])
    {
      pq = p;
    }
    else if(pqBag == tree[q])
    {
      pq = q;
    }
    else
    {
      pq = boost::add_vertex(pqBag, tree);
    }

    //Extend p-q edge with path from definition i.e. p - extension - q
    //Build extension from the middle pq vertex and outwards
    auto buildExtension =
      [&tree, pq]
      (const std::vector<unsigned int> & addition, unsigned int finalVertex)
    {
      auto previous = pq;
      for(auto v : addition)
      {
        auto newBag = tree[previous];
        newBag.insert(v);

        //Ignore last addition corresponding to 'finalVertex'
        if(newBag == tree[finalVertex])
        {
          break;
        }

        auto ext = boost::add_vertex(newBag, tree);
        boost::add_edge(ext, previous, tree);
        previous = ext;
      }

      if(previous != finalVertex)
      {
        boost::add_edge(previous, finalVertex, tree);
      }
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
      return boost::num_vertices(m_tree) == 1 ||
             boost::out_degree(vertex, m_tree) > 0;
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
 * Process the tree to satisfy critieria 3 of nice tree decompositions.
 *
 * Ensures that:
 * (1) all leafs have singleton bags
 * (2) the root node has an empty bag
 * 
 * @param inputTree Tree to be processed.
 * @param root Root node of the tree.
 * @return Processed tree satisfying critieria 3.
 */
static std::pair<tree_decomp_t, unsigned int> convertToCriteria3(
  const tree_decomp_t & inputTree,
  const unsigned int root)
{
  tree_decomp_t out(inputTree);

  //Extend vertex with a path containing decreasing bag sizes
  auto extendWithPath = [&out](unsigned int v)
  {
    unsigned int previous = v;
    while(out[previous].size() > 1)
    {
      auto bag = out[previous];
      bag.erase(bag.begin());
      auto next = boost::add_vertex(bag, out);
      boost::add_edge(previous, next, out);
      previous = next;
    }

    return previous;
  };

  //Generic BFS traversal calling the supplied callback
  auto traverse =
    [&inputTree, &out, root]
    (std::function<void(unsigned int, const std::list<unsigned int> &)> f)
  {
    std::queue<unsigned int> work;
    std::vector<bool> visited(boost::num_vertices(inputTree), false);

    work.push(root);
    visited.at(root) = true;
    while(!work.empty())
    {
      auto p = work.front(); work.pop();
      auto children = getChildren(out, p, visited);

      //Execute callback with current vertex and children list
      f(p, children);

      while(!children.empty())
      {
        work.push(children.front());
        visited.at(children.front()) = true;
        children.pop_front();
      }
    }
  };

  //Update all leaf nodes
  traverse(
    [&out, &extendWithPath]
    (unsigned int v, const std::list<unsigned int> & children)
  {
    //Process all leafs
    if(children.empty())
    {
      extendWithPath(v);
    }
  });

  //Update the root node
  unsigned int newRoot;
  traverse(
    [&out, &extendWithPath, &newRoot, root]
    (unsigned int v, const std::list<unsigned int> & children)
  {
    (void)children;
    //Process root node
    if(v == root)
    {
      newRoot = extendWithPath(v);

      if(!out[newRoot].empty())
      {
        auto previous = newRoot;
        newRoot = boost::add_vertex(std::set<unsigned int>(), out);;
        boost::add_edge(previous, newRoot, out);
      }
    }
  });

  assert(out[newRoot].empty());
  return std::make_pair(out, newRoot);
}

/**
 * Converts an existing valid tree decomposition into a nice tree decomposition.
 * A nice tree decomposition satisifies (Bodlaender and Kloks):
 * (1) Every node of the tree has at most two children
 * (2) if a node i has two children j; h then X_i = X_j = X_h,
 * (3) if a node i has one child, then either |X_i|=|X_j| + 1
 *     and X_j \subset X_i or |X_i|=|X_j| − 1 and X_i \subset X_j.
 *
 * @param inputTree Tree to be processed.
 * @return The resulting tree and a root node.
 */
std::pair<tree_decomp_t, unsigned int> convertToNiceDecomposition(const tree_decomp_t & inputTree)
{
  auto tree = convertToBinaryTree(inputTree);
  auto niceDecomp = convertToCriteria3(
                      convertToNiceCriteria2(
                        convertToNiceCriteria1(tree.first, tree.second),
                        tree.second),
                      tree.second);
  return niceDecomp;
}

} }
