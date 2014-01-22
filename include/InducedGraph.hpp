#pragma once

#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <set>

#include "TreeDecomposition.hpp"

namespace Count { namespace InducedGraph {

/**
 * Class implementing an edge predicate only keeping edges
 * (v_1, v_2) of which both v_1 and v_ 2 are in the provided set.
 */
template<typename T>
struct EdgePredicate
{
  EdgePredicate() = default;
  EdgePredicate(const T & g, const std::set<unsigned int> & f) : m_graph(g), m_filter(f) {}

  template<typename E>
  bool operator()(const E & e) const
  {
    return(m_filter.count(boost::source(e, m_graph)) > 0 && 
           m_filter.count(boost::target(e, m_graph)) > 0   );
  }

private:
 T m_graph;
 std::set<unsigned int> m_filter;
};

/**
 * Class implementing a vertex predicate only keeping vertices
 * that are in the provided set.
 */
template<typename T>
struct VertexPredicate
{
  VertexPredicate() = default;
  VertexPredicate(const T & g, const std::set<unsigned int> & f) : m_graph(g), m_filter(f) {}

  template<typename V>
  bool operator()(const V & v) const
  {
    return(m_filter.count(v) > 0);
  }

private:
 T m_graph;
 std::set<unsigned int> m_filter;
};

/**
 * Build an induced subgraph consisting of vertices in the given set.
 *
 * @param inputGraph Graph to be transformed.
 * @param acceptedVertices Set of vertices to be kept.
 * @return Induced subgraph.
 */
template<typename T> T buildInducedSubGraph(
  const T & inputGraph,
  const std::set<unsigned int> & acceptedVertices)
{
  EdgePredicate<T> edgeFilter(inputGraph, acceptedVertices);
  VertexPredicate<T> vertexFilter(inputGraph, acceptedVertices);

  boost::filtered_graph<T,
                        decltype(edgeFilter),
                        decltype(vertexFilter)>
                        filtered(inputGraph, edgeFilter, vertexFilter);

  T subGraph;
  boost::copy_graph(filtered, subGraph);
  return subGraph;
}


/**
 * Boost vertex copier ignoring any vertex properties.
 */
struct vertex_copier_ignore_prop {
  vertex_copier_ignore_prop() = default;

  template <typename Vertex1, typename Vertex2>
  void operator()(const Vertex1 & v1, Vertex2 & v2) const {
    (void)v1;
    (void)v2;
    return;
  }
};

/**
 * Build a tagged induced subgraph consisting of vertices in the given set.
 * The resulting vector is a map from old vertices onto new ones, or -1.
 *
 * @param inputGraph Graph to be transformed.
 * @param acceptedVertices Set of vertices to be kept.
 * @return Induced subgraph and the vertex map.
 */
template<typename T> std::pair<T, std::vector<int>> buildInducedSubGraphTagged(
  const T & inputGraph,
  const std::set<unsigned int> & acceptedVertices)
{
  typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS,
                              unsigned int>
          tagged_graph_t;

  auto const ignore_vertex_prop = boost::vertex_copy(vertex_copier_ignore_prop());

  //build tagged representation
  tagged_graph_t tagged;
  boost::copy_graph(inputGraph, tagged, ignore_vertex_prop);
  for(auto v = 0u; v < boost::num_vertices(tagged); ++v)
  {
    tagged[v] = v;
  }

  //apply filter on tagged graph
  EdgePredicate<tagged_graph_t>   edgeFilter(tagged, acceptedVertices);
  VertexPredicate<tagged_graph_t> vertexFilter(tagged, acceptedVertices);

  boost::filtered_graph<tagged_graph_t,
                        decltype(edgeFilter),
                        decltype(vertexFilter)>
                        filtered(tagged, edgeFilter, vertexFilter);

  //build map from old vertices onto new (or -1 if filtered)
  tagged_graph_t filteredCopy;
  boost::copy_graph(filtered, filteredCopy);

  std::vector<int> usage(boost::num_vertices(inputGraph), -1);
  for(decltype(filtered)::vertex_descriptor v = 0; v < boost::num_vertices(filteredCopy); ++v)
  {
    usage.at(filteredCopy[v]) = v;
  }

  //extract the filtered graph
  T subGraph;
  boost::copy_graph(filtered, subGraph, ignore_vertex_prop);

  return std::make_pair(subGraph, usage);
}

} }
