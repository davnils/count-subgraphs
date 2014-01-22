#include <set>
#include <tuple>

#include "Naive.hpp"

namespace Count { namespace Naive {

using namespace Utils;

/**
 *
 */
static unsigned long recurse(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host,
  const bool injective,
  std::set<unsigned int> patternSet,
  std::set<unsigned int> hostSet,
  std::map<unsigned int, unsigned int> & homo
  )
{
  if(patternSet.empty())
  {
    return 1;
  }

  auto const source = *patternSet.begin();
  patternSet.erase(patternSet.begin());
  unsigned long sum = 0;

  auto hostIt = hostSet.begin();
  while(hostIt != hostSet.end())
  {
    auto target = *hostIt;

    if(injective)
    {
      hostIt = hostSet.erase(hostIt);
    }
    else
    {
      ++hostIt;
    }

    //Check if valid extension of homomorphism
    bool isValid = true;
    for(auto it = boost::out_edges(source, pattern); it.first != it.second; ++it.first)
    {
      auto opposite = boost::target(*it.first, pattern);
      if(patternSet.count(opposite) == 0 &&
         !boost::edge(target, homo[opposite], host).second)
      {
        isValid = false;
        break;
      }
    }

    //Recurse on all valid extensions
    if(isValid)
    {
      homo[source] = target;
      sum += recurse(pattern, host, injective, patternSet, hostSet, homo);
      homo.erase(source);
    }

    if(injective)
    {
      hostSet.insert(hostIt, target);
    }
  }

  return sum;
}

/**
 *
 */
static unsigned long recursiveWrapper(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host,
  const bool injective)
{
  auto patternIt = boost::vertices(pattern);
  auto hostIt = boost::vertices(host);
  std::set<unsigned int> patternSet(patternIt.first, patternIt.second);
  std::set<unsigned int> hostSet(hostIt.first, hostIt.second);
  std::map<unsigned int, unsigned int> emptyHomo;

  return recurse(pattern, host, injective, patternSet, hostSet, emptyHomo);
}

/**
 *
 */
unsigned long countHomomorphisms(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host)
{
  return recursiveWrapper(pattern, host, false);
}

/**
 *
 */
unsigned long countInjective(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host)
{
  return recursiveWrapper(pattern, host, true);
}

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

std::ostream & operator<<(std::ostream & os, const std::vector<unsigned int> & set)
{
  os << " {";
  for(auto s : set)
  {
    os << (char)('a' + s) << " ";
  }
  os << "} ";

  return os;
};

std::ostream & operator<<(std::ostream & os, const std::map<unsigned int, unsigned int> & map)
{
  os << " {";
  for(auto e : map)
  {
    os << (char)('a' + e.first) << " -> " << (char)('a' + e.second) << ", ";
  }
  os << "} ";

  return os;
};

/**
 *
 */
unsigned long evaluateDisjointTriple(
  const std::pair<inj_homo_t, partition_triple_t> & homo,
  const Tree::undirected_graph_t & host,
  const PatternPartition & partition)
{
  //sum all products over all pairwise disjoint sets A,B,C of maximum size
  unsigned long result = 0;
  auto hostSet = std::set<unsigned int>(boost::vertices(host).first,
                                        boost::vertices(host).second);

  //remove vertices present in the homomorphism
  for(auto entry : homo.first)
  {
    hostSet.erase(entry.second);
  }

  auto counts = homo.second;
  auto outerEnum =
    [&hostSet, &counts, &result, &partition, &homo /* TODO: REMOVE*/]
    (const std::vector<bool> & markers)
  {
    auto vertexVec = extractSubset(markers, std::vector<unsigned int>(std::begin(hostSet), std::end(hostSet)));
    auto lVertices = std::set<unsigned int>(std::begin(vertexVec), std::end(vertexVec));

    //extract V(H) \ lVertices
    std::vector<unsigned int> mVertexRange;
    for(auto cand : hostSet)
    {
      if(lVertices.count(cand) == 0)
      {
        mVertexRange.push_back(cand);
      }
    }

    auto innerEnum =
      [&hostSet, &counts, &result, &lVertices, &mVertexRange, &homo, &partition]
      (const std::vector<bool> & markers)
    {
      //mVertices is extracted from V(H) \ L'
      auto vertexVec = extractSubset(markers, mVertexRange);
      auto mVertices = std::set<unsigned int>(std::begin(vertexVec), std::end(vertexVec));
      auto lmCombined = mergeSets(lVertices, mVertices);

      //extract range of vertices for R
      std::vector<unsigned int> rVertexRange;
      for(auto cand : hostSet)
      {
        if(lmCombined.count(cand) == 0)
        {
          rVertexRange.push_back(cand);
        }
      }

      auto innerMostEnum =
        [&hostSet, &counts, &result, &lVertices, &rVertexRange, &mVertices, &homo, &partition]
        (const std::vector<bool> & markers)
      {
        auto vertexVec = extractSubset(markers, rVertexRange);
        auto rVertices = std::set<unsigned int>(std::begin(vertexVec), std::end(vertexVec));

        if(counts.f.count(lVertices) == 0 ||
           counts.g.count(mVertices) == 0 ||
           counts.h.count(rVertices) == 0   )
         {
           return;
         }

        //Evaluate the triple
        result += counts.f.at(lVertices) *
                  counts.g.at(mVertices) *
                  counts.h.at(rVertices);
      };

      applyOnSubsets(rVertexRange.size(), partition.R.size(), innerMostEnum);
    };

    //choose |M| elements from V(H) \ L'
    applyOnSubsets(mVertexRange.size(), partition.M.size(), innerEnum);
  };

  //choose |L| elements from V(H)
  applyOnSubsets(hostSet.size(), partition.L.size(), outerEnum);
  return result;
}

/**
 *
 */
unsigned long evaluateDisjointTriples(
  const std::map<inj_homo_t, partition_triple_t> & counts,
  const Tree::undirected_graph_t & host,
  const PatternPartition & partition)
{
  unsigned long sum = 0;

  //for each injective homomorphism phi, evaluate the sum of products
  for(auto const & homo : counts)
  {
    auto term = evaluateDisjointTriple(homo, host, partition);
    sum += term;
  }

  return sum;
}

} }
