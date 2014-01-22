#include <algorithm>

#include "BasicEquation.hpp"

namespace Count { namespace Basic {

using namespace Utils;

/**
 * Calculate the intersection transform as defined in the counting paper.
 *
 * @param n Size of the universe.
 * @param z Z set.
 * @param q Size of the subset candidates.
 * @param t Current intersection size.
 * @param f Callback to be executed with an accepted subset.
 * @return Total sum.
 */
long calculateIntersectionTransform(
  const unsigned int n,
  const std::set<unsigned int> & z,
  const unsigned int q,
  const unsigned int t,
  std::function<long(const std::set<unsigned int> &)> f
  )
{
  long sum = 0;

  auto aggregate = [&sum, n, &z, f, t](const std::vector<bool> & markers)
  {
    auto const cand = extractSubset(markers);
    auto const candSet = std::set<unsigned int>(std::begin(cand), std::end(cand));

    assert(t <= n);

    if(buildSetIntersection<unsigned int>(z, candSet).size() == t)
    {
      sum += f(candSet);
    }
  };

  applyOnSubsets(n, q, aggregate);
  return sum;
}

/**
 * Calculate the parity transform as defined in the counting paper,
 * by calling the intersection transform.
 *
 * @param n Size of the universe.
 * @param z Z set.
 * @param q Size of the subset candidates.
 * @param s Range of intersection transform.
 * @param p Parity being evaluated.
 * @param f Callback to be executed with an accepted subset.
 * @return Total sum.
 */
long calculateParityTransform(
  const unsigned int n,
  const std::set<unsigned int> & z,
  const unsigned int q,
  const unsigned int s,
  const unsigned int p,
  std::function<long(const std::set<unsigned int> &)> f
  )
{
  long sum = 0;

  for(unsigned int t = 0; t <= s; ++t)
  {
    if(t % 2 == p)
    {
      sum += calculateIntersectionTransform(n, z, q, t, f);
    }
  }

  return sum;
}


/**
 * Calculate all tuples values as described in the counting paper.
 *
 * @param n Size of the universe.
 * @param q Size of the subset candidates.
 * @return Table of values indexed by i and s.
 */
std::vector<std::vector<unsigned long long>> calculateAllTuples(
  const unsigned int n,
  const unsigned int q
  )
{
  auto iMax = std::floor(3*q/2.0f);
  std::vector<std::vector<unsigned long long>> output
    (iMax+1, std::vector<unsigned long long>(2*iMax+1));

  for(auto i = 0u; i <= iMax; ++i)
  {
    for(auto s = 0u; s <= 2*iMax - i; ++s)
    {
      if(i == 0 && s == 0)
      {
        output.at(i).at(s) = 1;
      }
      else if(i < s)
      {
        output.at(i).at(s) = 0;
      }
      else if(i >= 1 && s == 0)
      {
        output.at(i).at(s) = output.at(i - 1).at(1);
      }
      else
      {
        output.at(i).at(s) =
          (n - s + 1) * output.at(i - 1).at(s - 1) +
          (s + 1)     * output.at(i - 1).at(s + 1);
      }
    }
  }

  return output;
}

/**
 * Calculate all T_p values as described in the counting paper.
 *
 * @param triple Triple containing partial counts.
 * @param q Size of the subset candidates.
 * @param n Size of the universe.
 * @param gamma Gamma parameter value being used.
 * @return Both parity vectors.
 */
std::pair<parity_vec_t, parity_vec_t> buildTpValues(
  const partition_triple_t & triple,
  const unsigned int q,
  const unsigned int n,
  const double gamma
  )
{
  unsigned int bound = std::floor((3/2.0f - gamma)*q);

  auto defaultCons = std::make_pair(parity_vec_t(bound+1), parity_vec_t(bound+1));
  parity_t f(defaultCons), g(defaultCons), h(defaultCons);

  //evaluate all parity transforms (for all s and p)
  for(unsigned int p = 0; p <= 1; ++p)
  {
    for(unsigned int s = 0; s <= bound; ++s)
    {
      //enumerate all subsets Z
      auto evalSubset =
        [&f, &g, &h, p, s, q, n, &triple]
        (const std::vector<bool> & markers)
      {
        auto subVec = extractSubset(markers);
        auto subset = std::set<unsigned int>(std::begin(subVec), std::end(subVec));

        auto fWrapper = [&triple](const std::set<unsigned int> & index)
        {
          if(triple.f.count(index) == 0)
          {
            return 0UL;
          }
          return triple.f.at(index);
        };

        auto gWrapper = [&triple](const std::set<unsigned int> & index)
        {
          if(triple.g.count(index) == 0)
          {
            return 0UL;
          }
          return triple.g.at(index);
        };

        auto hWrapper = [&triple](const std::set<unsigned int> & index)
        {
          if(triple.h.count(index) == 0)
          {
            return 0UL;
          }
          return triple.h.at(index);
        };

        if(p == 0)
        {
          f.first.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, fWrapper);
          g.first.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, gWrapper);
          h.first.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, hWrapper);
        }
        else
        {
          f.second.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, fWrapper);
          g.second.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, gWrapper);
          h.second.at(s)[subset] = calculateParityTransform(n, subset, q, s, p, hWrapper);
        }
      };
      applyOnSubsets(n, s, evalSubset);
    }
  }

  //Evaluate all T_0(Z) and T_1(Z), for all s and Z
  parity_vec_t T0(bound+1), T1(bound+1);
  for(unsigned int s = 0; s <= bound; ++s)
  {
    //enumerate all subsets Z
    auto evalSubset =
      [&f, &g, &h, s, &T0, &T1]
      (const std::vector<bool> & markers)
    {
      auto z = extractSubset(markers);
      auto zSet = std::set<unsigned int>(std::begin(z), std::end(z));

      auto pi = [&zSet, s](const parity_t & lookup, unsigned int p)
      {
        if(p == 0)
        {
          return lookup.first.at(s).at(zSet);
        }
        else
        {
          return lookup.second.at(s).at(zSet);
        }
      };

      T0.at(s)[zSet] =
        pi(f, 0) * pi(g, 0) * pi(h, 0) + 
        pi(f, 1) * pi(g, 1) * pi(h, 0) + 
        pi(f, 1) * pi(g, 0) * pi(h, 1) + 
        pi(f, 0) * pi(g, 1) * pi(h, 1);

      T1.at(s)[zSet] =
        pi(f, 1) * pi(g, 1) * pi(h, 1) + 
        pi(f, 0) * pi(g, 0) * pi(h, 1) + 
        pi(f, 0) * pi(g, 1) * pi(h, 0) + 
        pi(f, 1) * pi(g, 0) * pi(h, 0);
    };
    applyOnSubsets(n, s, evalSubset);
  }

  return std::make_pair(T0, T1);
}

/**
 * Build system of basic equations (y vector).
 *
 * @param triple Triple containing partial counts.
 * @param q Size of the subset candidates.
 * @param n Size of the universe.
 * @aram gamma Gamma parameter value being used.
 * @return Vector of y values.
 */
std::vector<long long> buildSystem(
  const partition_triple_t & triple,
  const unsigned int q,
  const unsigned int n,
  const double gamma
  )
{
  //evaluate all T_0 and T_1 values
  auto Tp = buildTpValues(triple, q, n, gamma);

  //evaluate all L_n
  auto L_n = calculateAllTuples(n, q);

  //evaluate all y_i
  std::vector<long long> yVec;
  unsigned int bound = std::floor((3/2.0f - gamma)*q);
  for(auto i = 0u; i <= bound; ++i)
  {
    long long total = 0;
    for(unsigned int s = 0; s <= i; ++s)
    {
      //enumerate all subsets Z
      long long tSum = 0;
      auto evalSubset =
        [s, &Tp, &tSum]
        (const std::vector<bool> & markers)
      {
        auto z = extractSubset(markers);
        auto zSet = std::set<unsigned int>(std::begin(z), std::end(z));
        tSum += Tp.first.at(s).at(zSet) - Tp.second.at(s).at(zSet);
      };
      applyOnSubsets(n, s, evalSubset);

      total += (long long)L_n.at(i).at(s) * tSum / (long long)binomCoefficient(n, s);
    }

    yVec.push_back(total);
  }

  return(yVec);
}

} }
