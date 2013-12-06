#include <boost/graph/random.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/random.hpp>

#include "TestUtils.hpp"
#include "TreeDecomposition.hpp"

namespace count { namespace test {

undirected_graph_t generateConnectedGraph(boost::random::mt19937 & gen)
{
  auto v = (rand() % 100)+1;
  auto e = (rand() % 100)+1;
  if(e > (v*(v-1))/2)
  {
    v = e;
  }

  //Generate a random graph
  undirected_graph_t graph;
  boost::generate_random_graph(graph, v, e, gen, false, false);

  std::vector<unsigned int> indexToComponent(boost::num_vertices(graph));
  auto components = boost::connected_components(graph, &indexToComponent[0]);

  //Connect all components, map each component to a vertex and add edges
  std::vector<unsigned int> componentMap(components);
  for(unsigned int i = 0; i < v; i++)
  {
    componentMap.at(indexToComponent.at(i)) = i;
  }

  for(unsigned int i = 0; i < components - 1; i++)
  {
    boost::add_edge(componentMap.at(i), componentMap.at(i + 1), graph);
  }

  return graph;
}

} }
