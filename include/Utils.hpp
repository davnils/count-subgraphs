#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <functional>
#include <iostream>
#include <map>
#include <vector>

#include "TreeDecomposition.hpp"

namespace Count { namespace Utils {

const unsigned int MAX_PATHWIDTH = 10;

typedef std::map<unsigned int, unsigned int> inj_homo_t;

struct partition_triple_t
{
  std::map<std::set<unsigned int>, unsigned long> f, g, h;
};

typedef boost::numeric::ublas::matrix<long long> matrix_t;

class PatternPartition
{
public:
  PatternPartition() = delete;
  explicit PatternPartition(
    const Tree::undirected_graph_t &
    );
  PatternPartition(
    const Tree::undirected_graph_t &,
    std::function<bool(PatternPartition*)> f
    );

  std::set<unsigned int> L, S, M, T, R;

  bool partitionExists() const { return found; }
  unsigned int getPathWidth() const { return pathWidth; }

private:
  bool found;
  unsigned int pathWidth;
  std::function<bool(PatternPartition*)> accept;

  void initialize(const Tree::undirected_graph_t &);

  bool recurse(
    const Tree::undirected_graph_t &, 
    const unsigned int,
    const unsigned int
    );
};

std::vector<unsigned int> extractSubset(const std::vector<bool> &);

std::vector<unsigned int> extractSubset(
  const std::vector<bool> &,
  const std::vector<unsigned int> &
  );

void applyOnSubsets(
  const unsigned int,
  const unsigned int,
  std::function<void(const std::vector<bool> &)>
  );

std::set<unsigned int> mergeSets(
  const std::set<unsigned int> &,
  const std::set<unsigned int> &
  );

/**
 * Extract the keys of a map.
 *
 * @param map The map.
 * @return Set of keys.
 */
template<typename K, typename V> std::set<K> extractKeys(const std::map<K, V> & map)
{
  std::set<K> keys;
  for(auto entry : map)
  {
    keys.insert(entry.first);
  }

  return keys;
}

/**
 * Calculate the bionomial coefficient.
 *
 * @param n n in "n choose k".
 * @param k k in "n choose k".
 * @return The result.
 */
template<typename T> T binomCoefficient(
  const T & n,
  const T & k
  )
{
  auto factorial = [](T a)
  {
    T prod = 1;

    while(a != 0)
    {
      prod *= a;
      --a;
    }

    return prod;
  };

  return factorial(n) / (factorial(n - k) * factorial(k));
}

/**
 * Take the intersection of two sets.
 *
 * @param a The first set.
 * @param b The second set.
 * @return a `intersection` b
 */
template<typename T> std::set<T> buildSetIntersection(
  const std::set<T> & a,
  const std::set<T> & b
  )
{
  std::vector<T> output(std::max(a.size(), b.size()));

  auto intersectEnd = std::set_intersection(std::begin(a),
                                            std::end(a),
                                            std::begin(b),
                                            std::end(b),
                                            std::begin(output));

  return std::set<T>(std::begin(output), intersectEnd);
}

/**
 * Take the difference of two sets.
 *
 * @param a The first set.
 * @param b The second set.
 * @return a \ b
 */
template<typename T> std::set<T> buildSetDifference(
  const std::set<T> & a,
  const std::set<T> & b
  )
{
  std::vector<T> output(std::max(a.size(), b.size()));

  auto differenceEnd = std::set_difference(std::begin(a),
                                           std::end(a),
                                           std::begin(b),
                                           std::end(b),
                                           std::begin(output));

  return std::set<T>(std::begin(output), differenceEnd);
}

/**
 * Take the symmetric difference of three sets.
 *
 * @param a The first set.
 * @param b The second set.
 * @param c The third set.
 * @return a `symdiff` b `symdiff` c
 */
template<typename T> std::set<T> calculateSymDiff(
  const std::set<T> & a,
  const std::set<T> & b,
  const std::set<T> & c
  )
{
  std::set<T> out;

  auto candidates = mergeSets(a, mergeSets(b, c));
  for(auto const & s : candidates)
  {
    if(((a.count(s) + b.count(s) + c.count(s)) % 2) == 1)
    {
      out.insert(s);
    }
  }

  return out;
}

/**
 * Take the symmetric difference of two sets.
 *
 * @param a The first set.
 * @param b The second set.
 * @return a `symdiff` b
 */
template<typename T> std::set<T> calculateSymDiff(
  const std::set<T> & a,
  const std::set<T> & b
  )
{
  return calculateSymDiff<T>(a, b, std::set<T>());
}

} }
