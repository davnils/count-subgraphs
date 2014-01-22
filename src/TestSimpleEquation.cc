#include <boost/numeric/ublas/matrix.hpp>
#include <boost/random.hpp>
#include <iostream>

#include "BasicEquation.hpp"
#include "Count.hpp"
#include "Naive.hpp"
#include "SimpleEquation.hpp"
#include "TestUtils.hpp"
#include "Utils.hpp"

namespace Count { namespace Simple { namespace Test {

using namespace Utils;
using namespace Count::Test;

std::ostream & operator<<(std::ostream & os, const std::set<unsigned int> & set)
{
  os << " {";
  for(auto s : set)
  {
    os << (char)('a' + s) << " ";
  }
  os << "} ";

  return os;
};

std::ostream & operator<<(std::ostream & os, const std::vector<int> & vec)
{
  os << " {";
  for(auto s : vec)
  {
    os << s << " ";
  }
  os << "} ";

  return os;
};

std::ostream & operator<<(std::ostream & os, const std::map<unsigned int, unsigned int> & map)
{
  os << " {";
  for(auto s : map)
  {
    os << (char)('a' + s.first)  << ":" << (char)('a' + s.second) << " ";
  }
  os << "} ";

  return os;
};

static unsigned long long naiveSymDiffProduct(
  const partition_triple_t & triple,
  const std::set<unsigned int> & dSet,
  const unsigned int n,
  const unsigned int q
  )
{
  unsigned long long sum = 0;

  auto generateA =
    [n, q, &triple, &dSet, &sum]
    (const std::vector<bool> & markers)
  {
    auto vec = extractSubset(markers);
    auto aSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

    auto generateA =
      [n, q, &triple, &dSet, &aSet, &sum]
      (const std::vector<bool> & markers)
    {
      auto vec = extractSubset(markers);
      auto bSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

      if(calculateSymDiff<unsigned int>(aSet, bSet) == dSet)
      {
        sum += triple.f.at(aSet) * triple.g.at(bSet);
      }
    };
    applyOnSubsets(n, q, generateA);
  };
  applyOnSubsets(n, q, generateA);

  return sum;
}

/**
 *
 */
static std::map<std::set<unsigned int>, unsigned long> generateSubsets(
  const unsigned int n,
  const unsigned int q,
  const unsigned int maxVal,
  boost::random::mt19937 & gen
  )
{
  std::map<std::set<unsigned int>, unsigned long> subsets;
  boost::random::uniform_int_distribution<> countDist(0, maxVal);

  //enumate all q-sized subsets taken of universe n
  auto enumerate =
    [&countDist, &gen, &subsets]
    (const std::vector<bool> & markers)
  {
    auto vec = extractSubset(markers);
    auto subset = std::set<unsigned int>(std::begin(vec), std::end(vec));

    //store random counts
    subsets[subset] = countDist(gen);
  };
  applyOnSubsets(n, q, enumerate);

  return subsets;
}

/**
 *
 */
void testSymDiffProduct(std::ostream & os)
{
  const unsigned int SEED = 0xf00d;
  boost::random::mt19937 gen(SEED);

  auto testCase = [&gen, &os]()
  {
    boost::random::uniform_int_distribution<> qDist(1, 3);
    unsigned int q = qDist(gen), n = 3*q;

    partition_triple_t triple;
    triple.f = generateSubsets(n, q, 5, gen);
    triple.g = generateSubsets(n, q, 5, gen);

    boost::random::uniform_int_distribution<> lDist(1, q);
    auto l = lDist(gen);
    l += (l % 2);

    std::cerr << "n=" << n << ", q=" << q << ", l=" << l << std::endl;

    auto diffSample = calculateDiffTransform(triple, q, n, l, l);

    for(auto & entry : diffSample)
    {
      auto ref = naiveSymDiffProduct(triple, entry.first, n, q);
      assert(ref == entry.second);
    }
  };

  runTests(testCase, "SymDiffProduct", os);
}

/**
 * Implementation of definition 2.5.
 */
long long naiveSolveSimpleIndeterminate1(
  const partition_triple_t & triple,
  const unsigned int n,
  const unsigned int q,
  const unsigned int j
  )
{
  long long result = 0;

  //enumerate all subsets A
  auto evalA =
    [n, q, j, &triple, &result]
    (const std::vector<bool> & markers)
  {
    auto vec = extractSubset(markers);
    auto aSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

    //enumerate all subsets B
    auto evalB =
      [n, q, j, &triple, &aSet, &result]
      (const std::vector<bool> & markers)
    {
      auto vec = extractSubset(markers);
      auto bSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

      //enumerate all subsets C
      auto evalC =
        [n, q, j, &triple, &aSet, bSet, &result]
        (const std::vector<bool> & markers)
      {
        auto vec = extractSubset(markers);
        auto cSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

        if(calculateSymDiff<unsigned int>(aSet, bSet, cSet).size() == j)
        {
          if(triple.f.count(aSet) != 0 &&
             triple.g.count(bSet) != 0 &&
             triple.h.count(cSet) != 0)
          {
            result += triple.f.at(aSet) *
                      triple.g.at(bSet) *
                      triple.h.at(cSet);
          }
        }
      };
      applyOnSubsets(n, q, evalC);
    };
    applyOnSubsets(n, q, evalB);
  };
  applyOnSubsets(n, q, evalA);

  return result;
}

/**
 * Implementation of definition 2.13.
 */
long long naiveSolveSimpleIndeterminate2(
  const partition_triple_t & triple,
  const unsigned int n,
  const unsigned int q,
  const unsigned int j
  )
{
  long long result = 0;

  //consider every l in the given range
  for(auto l = q - j; l <= q + j; ++l)
  {
    //enumerate all subsets 'D'
    auto evalD =
      [l, q, n, j, &triple, &result]
      (const std::vector<bool> & markers)
    {
      auto vec = extractSubset(markers);
      auto dSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

      //enumerate all subsets A and B
      long long sumAB = 0;
      auto evalA =
        [n, q, &triple, &dSet, &sumAB]
        (const std::vector<bool> & markers)
      {
        auto vec = extractSubset(markers);
        auto aSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

        auto evalB =
          [n, q, &triple, &dSet, aSet, &sumAB]
          (const std::vector<bool> & markers)
        {
          auto vec = extractSubset(markers);
          auto bSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

          if(buildSetDifference<unsigned int>(aSet, bSet) == dSet)
          {
            if(triple.f.count(aSet) != 0 && triple.g.count(bSet) != 0)
            {
              sumAB += triple.f.at(aSet) * triple.g.at(bSet);
            }
          }
        };
        applyOnSubsets(n, q, evalB);
      };
      applyOnSubsets(n, q, evalA);

      //enumerate all subsets C
      long long sumC = 0;
      auto evalC =
        [l, q, n, j, &triple, &dSet, &sumC]
        (const std::vector<bool> & markers)
      {
        auto vec = extractSubset(markers);
        auto cSet = std::set<unsigned int>(std::begin(vec), std::end(vec));

        assert((q + l - j) % 2 == 0);
        if(buildSetIntersection<unsigned int>(cSet, dSet).size() == (q + l - j)/2)
        {
          if(triple.h.count(cSet) != 0)
          {
            sumC += triple.h.at(cSet);
          }
        }
      };
      applyOnSubsets(n, q, evalC);

      result += sumAB * sumC;
    };
    applyOnSubsets(n, l, evalD);
  }

  return result;
}


/**
 *
 */
void testSimpleEquations(std::ostream & os)
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
        g1 = generateConnectedGraph(gen, 8);
      } while(boost::num_vertices(g1) < 3);

      do
      {
        g2 = generateConnectedGraph(gen, 8);
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

    auto injectiveRef = Naive::countInjective(g1,  g2);
    auto subResults = countSubgraphTriples(g1, g2, part);

    os << "Naive::countInjective()      : P->H " << injectiveRef << std::endl
       << "number of deltas generated: " << subResults.counts.size() << std::endl;

    auto const q = part.L.size();
    auto const n = boost::num_vertices(g2);

    const double gamma = 0.5f;
    auto const e = std::floor(3.0f*q/2.0f) + 1;
    auto const d = std::floor((3.0f/2.0f - gamma)*q) + 1;
    auto const simpleCount = e - d;

    //build vector containing all indeterminates (or rather, indices of)
    std::vector<unsigned int> indeterminates;
    for(unsigned int j = 0u; indeterminates.size() < simpleCount; ++j)
    {
      if(j % 2 == q % 2)
      {
        indeterminates.push_back(j);
      }
    }

    std::vector<unsigned int> indeterminatesBasic;
    for(unsigned int j = 0u; indeterminatesBasic.size() < e; ++j)
    {
      if(j % 2 == q % 2)
      {
        indeterminatesBasic.push_back(j);
      }
    }

    std::cerr << "simpleCount=" << simpleCount
              << ", total=" << e
              << std::endl;

    //disregard the case of insufficient amount of indeterminates
    if(simpleCount == 0)
    {
      return;
    }

    //solve for all suitable indeterminates, over all injective candidates
    int count = 1;
    for(auto const & triple : subResults.counts)
    {
      auto samples = solveIndeterminates(triple.second, q, n, simpleCount, gamma);

      //build basic-eq reference
      auto equations = Basic::buildSystem(triple.second, q, n, 0.0);
      boost::numeric::ublas::matrix<int> A(equations.second.size(), indeterminatesBasic.size()),
                                         y(equations.second.size(), 1);

      for(auto i = 0u; i < equations.second.size(); ++i)
      {
        for(auto j = 0u; j < indeterminatesBasic.size(); ++j)
        {
          auto var = indeterminatesBasic.at(j);
          auto coeff = Count::getCoefficient(n, i, var);
          A(i, j) = coeff;
          //os << coeff << " ";
        }

        y(i, 0) = equations.second.at(i);
        //os << " " << equations.second.at(i) << std::endl;
      }

      auto basicSolutions = Count::solveLinearSystem(A, y);
      //assert(basicSolutions.back() == Naive::evaluateDisjointTriple(triple, g2, part));

      //verify each indeterminate for this candidate (both through naive and basic-eq solvers)
      for(auto i = 0u; i < indeterminates.size(); ++i)
      {
        auto ref1 = naiveSolveSimpleIndeterminate1(triple.second, n, q, indeterminates.at(i));
        //auto ref2 = naiveSolveSimpleIndeterminate2(triple.second, n, q, indeterminates.at(i));
        auto sample = samples.at(i);

        /*std::cerr << "basic=" << basicSolutions.at(i)
                  << ", ref1=" << ref1
                  //<< ", ref2=" << ref2
                  << ", sample=" << sample
                  << std::endl;*/

        //assert(ref1 == ref2);
        assert(ref1 == basicSolutions.at(i));
        assert(sample == ref1);
      }

      count++;
    }
  };

  runTests(testCase, "SimpleEquations", os);
}

} } }
