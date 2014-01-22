#include <algorithm>
#include <boost/graph/graphviz.hpp>

#include "Utils.hpp"

namespace Count { namespace Utils {

/**
 * extract integers corresponding to the marked indices
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
 * extract integers from vector corresponding to the marked indices
 */
std::vector<unsigned int> extractSubset(
  const std::vector<bool> & markers,
  const std::vector<unsigned int> & from)
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
 * apply a function on all subsets of length 'length' with 'setCount' bits set
 */
void applyOnSubsets(
  const unsigned int length,
  const unsigned int setCount,
  std::function<void(const std::vector<bool> &)> f)
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

void visualizeGraph(
  std::ostream & os,
  const Tree::undirected_graph_t & graph)
{
  std::vector<std::string> labels =
    {"a", "b", "c", "d", "e", "f", "g",
     "h", "i", "j", "k", "l", "m", "n",
     "o", "p", "q", "r", "s", "t", "v",
     "w", "x", "y",  "z"};
  boost::write_graphviz(os, graph, boost::make_label_writer(&labels[0]));
  os << std::endl;
}

/**
 *
 */
std::set<unsigned int> mergeSets(
  const std::set<unsigned int> & s1,
  const std::set<unsigned int> & s2)
{
  std::set<unsigned int> copy(s1);
  copy.insert(std::begin(s2), std::end(s2));
  return copy;
};

/**
 *
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
 *
 */
PatternPartition::PatternPartition(
  const Tree::undirected_graph_t & graph,
  std::function<bool(PatternPartition*)> f
  ) : pathWidth(0), accept(f)
{
  initialize(graph);
}

/**
 *
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
 *
 * invariant: {L, S, M, T, R} always contains a valid (partial) partitioning.
 */
bool PatternPartition::recurse(
  const Tree::undirected_graph_t & graph, 
  const unsigned int pw,
  const unsigned int currentVertex)
{
  //Base case: all vertices have been assigned successfully
  if(currentVertex == boost::num_vertices(graph))
  {
    //TODO: Analyze if needed
    return accept(this);//S.size() && T.size();
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
