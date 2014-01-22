#include <algorithm>
#include <boost/graph/graphviz.hpp>
#include <boost/random.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>
#include <functional>
#include <iostream>

#include "Count.hpp"
#include "Naive.hpp"
#include "TestUtils.hpp"
#include "Utils.hpp"

namespace Count { namespace Test {

using namespace Utils;

/**
 * Verify the injective homomorphism counting procedure.
 * This is done by comparing against the naive implementation
 * and counting automorphisms.
 *
 * @param os Output stream to be written.
 */
void testCountInjective(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    auto g1 = generateConnectedGraph(gen, 3);
    decltype(g1) g2;

    do
    {
      g2 = generateConnectedGraph(gen, 5);
    } while(boost::num_vertices(g1) >= boost::num_vertices(g2));

    os << "|H|=" << boost::num_vertices(g2) << ", "
       << "|P|=" << boost::num_vertices(g1) << std::endl;

    unsigned long injectSum = 0;
    for(auto s : countInjective(g1, g2, std::set<unsigned int>(), inj_homo_t()))
    {
      injectSum += s.second;
    }

    auto injectRef = Naive::countInjective(g1,  g2);

    unsigned long autoCount = 0;
    for(auto s : countInjective(g1, g1, std::set<unsigned int>(), inj_homo_t()))
    {
      autoCount  += s.second;
    }

    auto autoRef = Naive::countInjective(g1, g1);

    os << "Count::countInjective(): P->H " << injectSum << std::endl
       << "Naive::countInjective(): P->H " << injectRef << std::endl
       << "Count::countInjective(): auto " << autoCount << std::endl
       << "Naive::countInjective(): auto " << autoRef   << std::endl;

    assert(injectSum == injectRef);
    assert(autoCount == autoRef);
    assert(injectSum % autoCount == 0);
  };

  runTests(testCase, "CountHomomorphisms", os);
}

/**
 * Verify that graph partitioning produces correct partitions, as
 * described in the counting paper.
 *
 * @param os Output stream to be written.
 */
void testPartitionGraph(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t graph;
    do
    {
      graph = generateConnectedGraph(gen, 13);
    } while(boost::num_vertices(graph) <= 2);

    os << "|G|=" << boost::num_vertices(graph) << std::endl;

    PatternPartition part(graph);
    if(!part.partitionExists())
    {
      std::cerr << "No partition found, skipping." << std::endl;
      return;
    }

    os
      << "|L| = " << part.L.size() << std::endl
      << "|S| = " << part.S.size() << std::endl
      << "|M| = " << part.M.size() << std::endl
      << "|T| = " << part.T.size() << std::endl
      << "|R| = " << part.R.size() << std::endl
      << "pw  = " << part.getPathWidth() << std::endl;

    auto lmrMax = std::ceil(boost::num_vertices(graph) / 3.0f);
    assert(part.L.size() <= lmrMax &&
           part.M.size() <= lmrMax &&
           part.R.size() <= lmrMax   );

    auto stMax = part.getPathWidth();
    assert(part.S.size() <= stMax &&
           part.T.size() <= stMax   );

    assert(part.L.size() +
           part.S.size() +
           part.M.size() +
           part.T.size() +
           part.R.size() == boost::num_vertices(graph));
  };

  runTests(testCase, "PartitionGraph", os);
}

/**
 * Verify counting injective homomorphisms using disjoint triples.
 * These are aggregated in a naive way and compared against a reference.
 *
 * @param os Output stream to be written.
 */
void testCountInjectiveTriples(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t g1, g2;
    do
    {
      g1 = generateConnectedGraph(gen, 6);
    } while(boost::num_vertices(g1) < 3);

    do
    {
      g2 = generateConnectedGraph(gen, 7);
    } while(boost::num_vertices(g1) > boost::num_vertices(g2));

    os << "|H| = " << boost::num_vertices(g2) << ", "
       << "|P| = " << boost::num_vertices(g1) << std::endl;

    auto predicate = 
      []
      (PatternPartition * p)
    {
        return p->S.size() > 0 &&  p->T.size() > 0 &&
               p->L.size() == p->M.size() && p->L.size() == p->R.size();
    };

    PatternPartition part(g1, predicate);

    if(!part.partitionExists())
    {
      std::cerr << "No partition found, skipping." << std::endl;
      return;
    }

    os
      << "|L| = " << part.L.size() << std::endl
      << "|S| = " << part.S.size() << std::endl
      << "|M| = " << part.M.size() << std::endl
      << "|T| = " << part.T.size() << std::endl
      << "|R| = " << part.R.size() << std::endl
      << "pw  = " << part.getPathWidth() << std::endl;

    auto injectRef = Naive::countInjective(g1,  g2);
    auto subResults = countSubgraphTriples(g1, g2, part);
    auto injectSample = Naive::evaluateDisjointTriples(subResults.counts, g2, part);

    os << "Count::countSubgraphTriples(): P->H " << injectSample << std::endl
       << "Naive::countInjective()      : P->H " << injectRef    << std::endl;

    assert(injectSample == injectRef);
  };

  runTests(testCase, "CountInjectiveTriples", os);
}

/**
 * Verify counting injective homomorphisms on a specific path.
 *
 * @param os Output stream to be written.
 */
void testCountInjectiveTriplesPath(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t g1(3), g2(3);

    boost::add_edge(0, 1, g1);
    boost::add_edge(1, 2, g1);

    boost::add_edge(0, 1, g2);
    boost::add_edge(1, 2, g2);

    os << "|H| = " << boost::num_vertices(g2) << ", "
       << "|P| = " << boost::num_vertices(g1) << std::endl;

    PatternPartition part(g1);

    os
      << "|L| = " << part.L.size() << std::endl
      << "|S| = " << part.S.size() << std::endl
      << "|M| = " << part.M.size() << std::endl
      << "|T| = " << part.T.size() << std::endl
      << "|R| = " << part.R.size() << std::endl
      << "pw  = " << part.getPathWidth() << std::endl;

    auto injectRef = Naive::countInjective(g1,  g2);
    auto subResults = countSubgraphTriples(g1, g2);
    auto injectSample = Naive::evaluateDisjointTriples(subResults.counts, g2, subResults.part);

    os << "Count::countSubgraphTriples() (ext) : P->H " << injectSample << std::endl
       << "Naive::countInjective()             : P->H " << injectRef    << std::endl;

    assert(injectSample == injectRef);
  };

  runTests(testCase, "CountInjectiveTriplesPath", os);
}

/**
 * Verify counting injective homomorphisms on k-paths.
 *
 * @param os Output stream to be written.
 */
void testCountInjectiveTriplesKPath(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    boost::random::uniform_int_distribution<> lengthDist(3, 9);
    unsigned int pathLength = lengthDist(gen);
    Tree::undirected_graph_t g1(pathLength), g2(pathLength);

    for(unsigned int v = 0; v < pathLength - 1; ++v)
    {
      boost::add_edge(v, v + 1, g1);
      boost::add_edge(v, v + 1, g2);
    }

    os << "path length: " << pathLength << std::endl;

    auto predicate = 
      []
      (PatternPartition * p)
    {
        return p->S.size() > 0 &&  p->T.size() > 0 &&
               p->L.size() == p->M.size() && p->L.size() == p->R.size();
    };

    PatternPartition part(g1, predicate);

    if(!part.partitionExists())
    {
      std::cerr << "No partition found, skipping." << std::endl;
      return;
    }

    os
      << "|L| = " << part.L.size() << std::endl
      << "|S| = " << part.S.size() << std::endl
      << "|M| = " << part.M.size() << std::endl
      << "|T| = " << part.T.size() << std::endl
      << "|R| = " << part.R.size() << std::endl
      << "pw  = " << part.getPathWidth() << std::endl;

    auto injectRef = Naive::countInjective(g1,  g2);
    auto subResults = countSubgraphTriples(g1, g2, part);
    auto injectSample = Naive::evaluateDisjointTriples(subResults.counts, g2, part);

    os << "Count::countInjective() (ext) : P->H " << injectSample << std::endl
       << "Naive::countInjective()       : P->H " << injectRef    << std::endl;

    assert(injectSample == injectRef);
  };

  runTests(testCase, "CountInjectiveTriplesKPath", os);
}

/**
 * Verify counting subgraph isomorphisms.
 * Compared against a reference by counting in a naive fashion.
 * Note that invalid partitions are simply skipped.
 *
 * @param os Output stream to be written.
 */
void testCountIsoSubgraphs(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t g1, g2;
    do
    {
      g1 = generateConnectedGraph(gen, 9);
    } while(boost::num_vertices(g1) < 3);

    do
    {
      g2 = generateConnectedGraph(gen, 9);
    } while(boost::num_vertices(g1) > boost::num_vertices(g2));

    os << "|H| = " << boost::num_vertices(g2) << ", "
       << "|P| = " << boost::num_vertices(g1) << std::endl;

    auto sample = countIsoSubgraphs(g1, g2);

    if(!sample.first)
    {
      std::cerr << "No partitioning found." << std::endl;
      return;
    }

    auto ref = Naive::countInjective(g1, g2) / Naive::countInjective(g1,  g1);

    os << "Count::countIsoSubgraphs()  :   P->H " << sample.second << std::endl
       << "Naive::countInjective()/auto:   P->H " << ref << std::endl;

    assert(sample.second == ref);
  };

  runTests(testCase, "CountIsoSubgraphs", os);
}

} }
