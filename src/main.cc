#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>

#include "TestTreeDecomposition.hpp"
#include "TreeDecomposition.hpp"

static void checkDiaz()
{
  Tree::undirected_graph_t graph(7);
  enum {A, B, C, D, E, F, G} edge_t;

  boost::add_edge(A, B, graph);
  boost::add_edge(A, G, graph);
  boost::add_edge(B, C, graph);
  boost::add_edge(B, G, graph);
  boost::add_edge(C, D, graph);
  boost::add_edge(C, E, graph);
  boost::add_edge(D, E, graph);
  boost::add_edge(E, F, graph);
  boost::add_edge(E, G, graph);
  boost::add_edge(F, G, graph);
  boost::add_edge(F, C, graph);

  std::vector<std::string> labels = {"a", "b", "c", "d", "e", "f", "g"};
  boost::write_graphviz(std::cout, graph, boost::make_label_writer(&labels[0]));
  std::cout << std::endl;

  auto result = Tree::buildTreeDecomposition(graph);
  Tree::visualizeDecomposition(std::cout, result);

  auto nice = Tree::convertToNiceDecomposition(result);
  Tree::visualizeDecomposition(std::cout, nice.first);
  std::cout << "Root vertex: " << nice.second << std::endl;
}

int main()
{
  //Tree::Test::testBinaryDecomposition(std::cout);
  Tree::Test::testNiceTreeDecomposition(std::cout);
  Tree::Test::testCountHomomorphisms(std::cout);
  return 0;
}
