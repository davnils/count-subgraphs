#include <boost/graph/random.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/random.hpp>

#include "TestUtils.hpp"
#include "TreeDecomposition.hpp"

namespace Count { namespace Test {

/**
 * Generate a connected graph of an approximate maximum size (number of vertices).
 *
 * @param gen Number generator to be used.
 * @param maxSize Approximate maximum number of vertices.
 */
Tree::undirected_graph_t generateConnectedGraph(
  boost::random::mt19937 & gen,
  const unsigned int maxSize
  )
{
  assert(maxSize >= 1);

  auto v = (rand() % maxSize)+1;
  auto e = (rand() % maxSize)+1;
  if(e > (v*(v-1))/2)
  {
    v = e + 1;
  }

  //Generate a random graph
  Tree::undirected_graph_t graph;
  boost::generate_random_graph(graph, v, e, gen, false, false);

  std::vector<unsigned int> indexToComponent(boost::num_vertices(graph));
  auto components = boost::connected_components(graph, &indexToComponent[0]);

  //Connect all components, map each component to a vertex and add edges
  std::vector<unsigned int> componentMap(components);
  for(unsigned int i = 0; i < v; ++i)
  {
    componentMap.at(indexToComponent.at(i)) = i;
  }

  for(unsigned int i = 0; i < components - 1; ++i)
  {
    boost::add_edge(componentMap.at(i), componentMap.at(i + 1), graph);
  }

  return graph;
}

/**
 * Execute a test case with the given description.
 *
 * @param testCase Test to be called.
 * @param title Title of the test.
 * @param os Output stream to be written.
 */
void runTests(
  std::function<void(void)> testCase,
  const std::string & title,
  std::ostream & os
  )
{
  const unsigned int TESTS = 100;

  for(unsigned int test = 0; test < TESTS; ++test)
  {
    os
      << "------------------NEW TEST---------------------\n"
      << "Running test (" << title << "): #" << test << "\n";

    testCase();

    os
      << "-----------------------------------------------\n"
      << std::endl;
  }
}


} }
