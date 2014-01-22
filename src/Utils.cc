#include <algorithm>
#include <boost/graph/graphviz.hpp>

#include "Utils.hpp"

namespace Count { namespace Utils {

/**
 * Extract integers corresponding to the marked indices.
 * @param markers True if the corresponding integer (position) should be included.
 * @return Vector containing all the integers.
 */
std::vector<unsigned int> extractSubset(const std::vector<bool> & markers)
{
  std::vector<unsigned int> items;
  for(auto i = 0u; i < markers.size(); ++i)
  {
    if(markers.at(i))
    {
      items.push_back(i);
    }
  }

  return items;
}

/**
 * Extract integers corresponding to the marked indices in an array.
 * @param markers True if the corresponding integer should be included.
 * @return Vector containing all the integers.
 */
std::vector<unsigned int> extractSubset(
  const std::vector<bool> & markers,
  const std::vector<unsigned int> & from
  )
{
  std::vector<unsigned int> items;
  for(auto i = 0u; i < markers.size(); ++i)
  {
    if(markers.at(i))
    {
      items.push_back(from.at(i));
    }
  }

  return items;
}

/**
 * Apply a function on all subsets with size decided by the parameters.
 *
 * @param length Total number of elements.
 * @param setCount Number of elements chosen in each subset.
 * @param f Function to be called on each subset.
 */
void applyOnSubsets(
  const unsigned int length,
  const unsigned int setCount,
  std::function<void(const std::vector<bool> &)> f
  )
{
  assert(setCount <= length && setCount >= 0);

  unsigned long long subset = (1ULL << setCount) - 1;

  auto gosperNext =
    []
    (unsigned long long & x)
  {
    auto u = x & -x;
    auto v = u + x;

    if(v == 0)
    {
      return false;
    }

    x = v + (((v^x) / u) >> 2);
    return true;
  };

  std::vector<bool> vec(length, false);
  do
  {
    for(auto i = 0u; i < length; ++i)
    {
      vec.at(i) = subset & (1ULL << i);
    }
    f(vec);
  } while(setCount && gosperNext(subset) && subset < (1ULL << length));
}

/**
 * Take the union of two sets.
 *
 * @param s1 The first set.
 * @param s2 The second set.
 * @return s1 `union` s2
 */
std::set<unsigned int> mergeSets(
  const std::set<unsigned int> & s1,
  const std::set<unsigned int> & s2
  )
{
  std::set<unsigned int> copy(s1);
  copy.insert(std::begin(s2), std::end(s2));
  return copy;
};

/**
 * Constructor accepting the graph to be partitioned.
 *
 * @param graph The input graph.
 */
PatternPartition::PatternPartition(
  const Tree::undirected_graph_t & graph
  ) : pathWidth(0)
{
  accept = 
    []
    (PatternPartition * p)
  {
      return p->S.size() && p->T.size();
  };

  initialize(graph);
}

/**
 * Constructor accepting the graph to be partitioned and a predicate
 * which defines when a candidate partition is accepted.
 *
 * @param graph The input graph.
 * @param f Predicate defining termination.
 */
PatternPartition::PatternPartition(
  const Tree::undirected_graph_t & graph,
  std::function<bool(PatternPartition*)> f
  ) : pathWidth(0), accept(f)
{
  initialize(graph);
}

/**
 * Calculate a 5-partition of the given graph.
 *
 * @param graph The input graph.
 */
void PatternPartition::initialize(const Tree::undirected_graph_t & graph)
{
  found = false;
  for(auto pw = 1u; pw <= MAX_PATHWIDTH; ++pw)
  {
    if((found = recurse(graph, pw, 0)))
    {
      pathWidth = pw;
      found = true;
      break;
    }
  }
}

/**
 * Recursive procedure deciding if there is some 5-partition
 * accepted by the predicate.
 *
 * @param graph The input graph.
 * @param pw The maximum pathwidth to be used.
 * @param currentVertex The vertex to be considered, or |G|.
 * @return True if the predicate accepted a partition.
 */
bool PatternPartition::recurse(
  const Tree::undirected_graph_t & graph, 
  const unsigned int pw,
  const unsigned int currentVertex
  )
{
  //Base case: all vertices have been assigned successfully
  if(currentVertex == boost::num_vertices(graph))
  {
    return accept(this);
  }

  //evaluate a candidate placement
  auto tryCand =
    [this, &graph, pw, currentVertex]
    (std::set<unsigned int> & vertexSet, const std::set<unsigned int> & allowed, unsigned maxSize)
  {
    //check if all edges associated with the vertex go to vertices in 'allowed'
    for(auto it = boost::adjacent_vertices(currentVertex, graph); it.first != it.second; ++it.first)
    {
      auto other = *it.first;
      if(other < currentVertex && allowed.count(other) == 0)
      {
        return false;
      }
    }

    //check if there's space available
    if(vertexSet.size() >= maxSize)
    {
      return false;
    }

    //place vertex into set
    vertexSet.insert(currentVertex);
    if(this->recurse(graph, pw, currentVertex + 1))
    {
      return true;
    }
    else
    {
      vertexSet.erase(currentVertex);
      return false;
    }
  };

  //Consider the five possible placements of the current vertex
  //Will only recurse with valid solutions

  auto lmrMaxCount = std::ceil(boost::num_vertices(graph)/3.0f);
  auto stMaxCount  = pw;
  if(tryCand(L, mergeSets(L, S), lmrMaxCount))
  {
    return true;
  }

  if(tryCand(S, mergeSets(L, mergeSets(S, M)), stMaxCount))
  {
    return true;
  }

  if(tryCand(M, mergeSets(S, mergeSets(M, T)), lmrMaxCount))
  {
    return true;
  }

  if(tryCand(T, mergeSets(M, mergeSets(T, R)), stMaxCount))
  {
    return true;
  }

  if(tryCand(R, mergeSets(T, R), lmrMaxCount))
  {
    return true;
  }

  return false;
}

} }
