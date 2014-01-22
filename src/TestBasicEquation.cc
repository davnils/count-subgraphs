#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <iostream>

#include "BasicEquation.hpp"
#include "Count.hpp"
#include "Naive.hpp"
#include "TestUtils.hpp"
#include "Utils.hpp"

namespace Count { namespace Basic { namespace Test {

using namespace Utils;
using namespace Count::Test;

/**
 * Verify the implementation of the L table recurrence, 
 * by comparing against a naive recursion.
 *
 * @param os Output stream to be written.
 */
void testLTable(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    boost::random::uniform_int_distribution<> qDist(1, 15);
    unsigned int q = qDist(gen), n = 3*q;

    std::function<unsigned long long(unsigned int, unsigned int)> naive;
    naive =
      [n, &naive]
      (unsigned int i, unsigned int s) -> unsigned long long
    {
      if(i == 0 && s == 0)
      {
        return 1;
      }
      else if(i < s)
      {
        return 0;
      }
      else if(i >= 1 && s == 0)
      {
        return naive(i - 1, 1);
      }
      else
      {
        return (n - s + 1) * naive(i - 1, s - 1) + 
               (s + 1)     * naive(i - 1, s + 1);
      }
    };

    auto tuples = calculateAllTuples(n, q);
    for(auto i = 0; i <= std::floor(3*q/2.0f); ++i)
    {
      for(auto s = 0; s <= i; ++s)
      {
        assert(naive(i, s) == tuples.at(i).at(s));
      }
    }
  };

  runTests(testCase, "LTable", os);
}

/**
 * Retrieve the subset indicated by the markers as a set.
 *
 * @param markers Markers indicating the integers to be saved.
 * @return Set of integers.
 */
static std::set<unsigned int> getSubset(const std::vector<bool> & markers)
{
  auto const vertices = extractSubset(markers);
  return std::set<unsigned int>(std::begin(vertices), std::end(vertices));
}

/**
 * Calculate T_p values (as described in the counting paper) by naively
 * enumerating all subsets.
 *
 * @param q Size of the subset candidates.
 * @param n Size of the universe.
 * @param p Parity value (0 or 1).
 * @param z Z set.
 * @param triple Triple containing partial counts.
 * @return Total sum.
 */
static long long naiveCalculateTp(
  const unsigned int q,
  const unsigned int n,
  const unsigned int p,
  const std::set<unsigned int> & z,
  const partition_triple_t & triple
  )
{
  long long sum = 0;

  //enumerate all q-sized subsets A,B,C
  auto enumerateA = 
    [q, n, p, &z, &triple, &sum]
    (const std::vector<bool> & markers)
  {
    auto const aVertices = getSubset(markers);

    auto enumerateB = 
      [q, n, p, &z, &triple, &sum, &aVertices]
      (const std::vector<bool> & markers)
    {
      auto const bVertices = getSubset(markers);

      auto enumerateC = 
        [q, n, p, &z, &triple, &sum, &aVertices, &bVertices]
        (const std::vector<bool> & markers)
      {
        auto const cVertices = getSubset(markers);

        //sum all entries satisfiying |symDiff(a,b,c) `intersect` z| ~= p
        auto diff = calculateSymDiff(aVertices, bVertices, cVertices);

        if(buildSetIntersection<unsigned int>(diff, z).size() % 2 == p)
        {
          if(triple.f.count(aVertices) != 0 &&
             triple.g.count(bVertices) != 0 &&
             triple.h.count(cVertices) != 0   )
          {
            sum += triple.f.at(aVertices) *
                   triple.g.at(bVertices) *
                   triple.h.at(cVertices);
          }
        }
      };
      applyOnSubsets(n, q, enumerateC);

    };

    applyOnSubsets(n, q, enumerateB);
  };

  applyOnSubsets(n, q, enumerateA);
  return sum;
}

/**
 * Naively calculate a specific y value.
 *
 * @param i The equation index, 0 to e.
 * @param q Size of the subset candidates.
 * @param n Size of the universe.
 * @param triple Triple containing partial counts.
 * @return Total sum.
 */
static long long naiveCalculateYi(
  const unsigned int i,
  const unsigned int q,
  const unsigned int n,
  const partition_triple_t & triple
  )
{
  long long sum = 0;

  auto nextNTuple =
    []
    (const unsigned int n, std::vector<unsigned int> & vec)
  {
    auto const length = vec.size();

    auto nextIndex = 0u;
    while(nextIndex < length &&
          vec.at(nextIndex) == n - 1)
    {
      vec.at(nextIndex) = 0;
      ++nextIndex;
    }

    if(nextIndex == length)
    {
      return false;
    }

    vec.at(nextIndex)++;
    return true;
  };

  auto tpSamples = buildTpValues(triple, q, n, 0.0);

  std::vector<unsigned int> tuple(i, 0);
  do
  {
    //keep all elements occuring an odd number of times
    std::set<unsigned int> z;
    for(auto cand : tuple)
    {
      if(std::count(std::begin(tuple), std::end(tuple), cand) % 2 == 1)
      {
        z.insert(cand);
      }
    }

    //add the diff of both T_p values
    auto t0 = naiveCalculateTp(q, n, 0, z, triple);
    auto t1 = naiveCalculateTp(q, n, 1, z, triple);

    if(t0 != tpSamples.first.at(z.size()).at(z))
    {
      assert(false);
    }

    if(t1 != tpSamples.second.at(z.size()).at(z))
    {
      assert(false);
    }

    sum += t0 - t1;
  } while(nextNTuple(n, tuple));

  return sum;
}


/**
 * Verify that basic equations solve for the expected injective count.
 *
 * @param os Output stream to be written.
 */
void testBasicEquations(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    Tree::undirected_graph_t g1, g2;

    auto predicate = 
      []
      (PatternPartition * p)
    {
        return p->S.size() > 0 &&  p->T.size() > 0 &&
               p->L.size() == p->M.size() && p->L.size() == p->R.size();
    };

    PatternPartition part(g1, predicate);

    do
    {
      do
      {
        g1 = generateConnectedGraph(gen, 6);
      } while(boost::num_vertices(g1) < 3);

      do
      {
        g2 = generateConnectedGraph(gen, 6);
      } while(boost::num_vertices(g1) > boost::num_vertices(g2));

      part = PatternPartition(g1, predicate);
    } while(!part.partitionExists());

    assert(part.L.size() == part.M.size() && part.L.size() == part.R.size());

    os << "|H| = " << boost::num_vertices(g2) << ", "
       << "|P| = " << boost::num_vertices(g1) << std::endl;

    os
      << "|L| = " << part.L.size() << std::endl
      << "|S| = " << part.S.size() << std::endl
      << "|M| = " << part.M.size() << std::endl
      << "|T| = " << part.T.size() << std::endl
      << "|R| = " << part.R.size() << std::endl
      << "pw  = " << part.getPathWidth() << std::endl;

    auto injectRef = Naive::countInjective(g1,  g2);
    auto subResults = countSubgraphTriples(g1, g2, part);

    os << "Naive::countInjective()      : P->H " << injectRef << std::endl
       << "number of deltas generated: " << subResults.counts.size() << std::endl;

    auto const q = part.L.size();

    //build vector containing all indeterminates (or rather, indices of)
    std::vector<unsigned int> indeterminates;
    for(unsigned int j = 0u; j <= 3*q; ++j)
    {
      if(j % 2 == (q % 2))
      {
        indeterminates.push_back(j);
      }
    }

    unsigned long long injectiveSum = 0;
    for(auto const & triple : subResults.counts)
    {
      //build the equation system
      auto equations = buildSystem(triple.second, q, boost::num_vertices(g2), 0.0);

      boost::numeric::ublas::matrix<int> A(equations.second.size(), indeterminates.size()),
                                         y(equations.second.size(), 1);

      for(auto i = 0u; i < equations.second.size(); ++i)
      {
        for(auto j = 0u; j < indeterminates.size(); ++j)
        {
          auto var = indeterminates.at(j);
          auto coeff = Count::getCoefficient(boost::num_vertices(g2), i, var);
          A(i, j) = coeff;
          os << coeff << " ";
        }

        y(i, 0) = equations.second.at(i);
        os << " " << equations.second.at(i);
        std::flush(os);

        auto eqRef = naiveCalculateYi(i, q, boost::num_vertices(g2), triple.second);
        os << " (naive: " << eqRef << ")";
        os << std::endl;
        assert(eqRef == equations.second.at(i));
      }

      //add up the contribution gathered by solving the system
      auto contrib = Count::solveLinearSystem(A, y).back();

      auto contribRef = Naive::evaluateDisjointTriple(triple, g2, part);

      assert(contrib >= 0);
      assert(static_cast<unsigned long>(contrib) == contribRef);

      injectiveSum += contrib;
    }

    os << "sum=" << injectiveSum << " and ref="<< injectRef << std::endl;
    assert(injectRef == injectiveSum);
  };

  runTests(testCase, "BasicEquations", os);
}

} } }
