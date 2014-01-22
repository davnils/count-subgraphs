#include <algorithm>
#include <boost/graph/graphviz.hpp>
#include <boost/random.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <functional>
#include <iostream>

#include "Homomorphism.hpp"
#include "Naive.hpp"
#include "TestUtils.hpp"
#include "TreeDecomposition.hpp"

namespace Count { namespace Homomorphism { namespace Test {

using namespace Count::Test;

/**
 * Test stingy orderings.
 * @param os Stream to be written with output.
 */
void testStingyOrdering(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    auto randomGraph = generateConnectedGraph(gen, 200);
    auto decomp = Tree::buildTreeDecomposition(randomGraph);
    auto niceTree = Tree::convertToNiceDecomposition(decomp);

    auto ordering = calculateStingyOrdering(niceTree.first, niceTree.second);

    //(1) all vertices are included exactly once
    assert(ordering.size() == boost::num_vertices(niceTree.first));
    std::vector<unsigned int> copy(std::begin(ordering), std::end(ordering));
    std::sort(std::begin(copy), std::end(copy));
    std::unique(std::begin(copy), std::end(copy));
    assert(copy.size() == ordering.size());

    //(2) the root vertex is in the last position
    assert(ordering.back() == niceTree.second);
 
    //(3) the first position is a child (out-degree 1, or 0 if singleton)
    assert(ordering.front() != niceTree.second);
    assert(boost::out_degree(ordering.front(), niceTree.first) <= 1);
  };

  runTests(testCase, "StingyOrdering", os);
}

/**
 * Test naive counting, i.e. exhaustive search.
 * Injective and non-injective counts are compared.
 * Also, a known number of automorphisms is verified.
 *
 * @param os Stream to be written with output.
 */
void testNaiveCounting(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    auto host = generateConnectedGraph(gen, 500);

    Tree::undirected_graph_t triangle(3);
    boost::add_edge(0, 1, triangle);
    boost::add_edge(1, 2, triangle);
    boost::add_edge(2, 0, triangle);

    const unsigned long triangleAutomorphisms = 6;

    assert(Naive::countInjective(triangle, triangle) == triangleAutomorphisms);

    auto injCount = Naive::countInjective(triangle, host);
    assert(injCount % triangleAutomorphisms == 0);

    auto homoCount = Naive::countHomomorphisms(triangle, host);
    assert(homoCount >= injCount);

    auto autoCount = Naive::countHomomorphisms(triangle, triangle);
    assert(autoCount == triangleAutomorphisms);

    os << "Triangle count: "
       << injCount
       << std::endl;
  };

  runTests(testCase, "NaiveInjectiveHomomorphism", os);
}

/**
 * Verify that counting homomorphisms gives some known value.
 *
 * @param os Stream to be written with output.
 */
void testSimpleHomomorphism(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t g1(2), g2(2);

    boost::add_edge(0, 1, g1);
    boost::add_edge(0, 1, g2);

    auto decomp = Tree::buildTreeDecomposition(g1);
    auto niceTree = Tree::convertToNiceDecomposition(decomp);

    map_t homo1 = { std::make_pair(0, 0) },
          homo2 = { std::make_pair(1, 0) };
    auto homoCount1 = countHomomorphisms(g1, niceTree.first, niceTree.second, g2, homo1);
    auto homoCount2 = countHomomorphisms(g1, niceTree.first, niceTree.second, g2, homo2);
    auto refCount = 1u;

    os << "Homomorphism::countHomomorphisms() : "
       << homoCount1 << " and " << homoCount2
       << std::endl;

    assert(homoCount1 == homoCount2 && homoCount1 == refCount);
  };

  runTests(testCase, "SimpleHomomorphism", os);
}

/**
 * Verify that the naive and proper implementation of
 * counting homomorphisms give the same results.
 *
 * @param os Stream to be written with output.
 */
void testCountHomomorphisms(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    auto g1 = generateConnectedGraph(gen, 7);
    decltype(g1) g2;

    do
    {
      g2 = generateConnectedGraph(gen, 12);
    } while(boost::num_vertices(g1) >= boost::num_vertices(g2));

    auto injective = Naive::countInjective(g1,  g2);
    os << "|H|=" << boost::num_vertices(g2) << ", "
       << "|P|=" << boost::num_vertices(g1)  << ", "
       << "count=" << injective << std::endl;
    assert(injective >= 0);

    auto decomp = Tree::buildTreeDecomposition(g1);
    auto niceTree = Tree::convertToNiceDecomposition(decomp);

    map_t emptyHomo;
    auto homoCount = countHomomorphisms(g1, niceTree.first, niceTree.second, g2, emptyHomo);
    auto refCount = Naive::countHomomorphisms(g1, g2);

    os << "Homomorphism::countHomomorphisms() : " << homoCount << std::endl;
    os << "Naive::countHomomorphisms()        : " << refCount << std::endl;

    assert(homoCount == refCount);
  };

  runTests(testCase, "CountHomomorphisms", os);
}

} } }
