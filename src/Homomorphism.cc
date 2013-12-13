#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <map>
#include <set>

#include "Homomorphism.hpp"
#include "TreeDecomposition.hpp"

namespace Homomorphism {

static std::list<unsigned int> calculateStingyOrdering(
  const Tree::tree_decomp_t & tree)
{
  std::list<unsigned int> stingy;
  //TODO
  return stingy;
}

/**
 *
 */
unsigned long long countHomomorphisms(
  const Tree::undirected_graph_t & pattern,
  const Tree::tree_decomp_t & decomp,
  const unsigned int root,
  const Tree::undirected_graph_t & host)
{
  typedef std::map<unsigned int, unsigned int> map_t;
  std::map<unsigned int, std::map<map_t, count_t>> counts;

  auto ordering = calculateStingyOrdering(decomp);
  std::vector<bool> visited(ordering.size(), false);

  auto getChilds = [&decomp, &visited](unsigned int v)
  {
    std::set<unsigned int> neighbours;
    for(auto it = boost::adjacent_vertices(v, decomp); it.first != it.second;
        ++it.first)
    {
      if(visited.at(v))
      {
        neighbours.insert(v);
      }
    }

    return neighbours;
  };

  auto nextMap = [&host](std::vector<unsigned int> & vec)
  {
    auto alphabetSize = boost::num_vertices(host);
    auto length = vec.size();

    auto nextIndex = 0;
    while(nextIndex < length &&
          vec.at(nextIndex) == alphabetSize - 1)
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

  auto getUniqueSetDiff = [&decomp](unsigned int first, unsigned int second,
                                    const std::set<unsigned int> & children)
  {
    std::vector<unsigned int> diff(1);
    std::set_difference(decomp[first].begin(), decomp[first].end(),
                        decomp[second].begin(), decomp[second].end(),
                        diff.begin());
    return diff.at(0);
  };

  //TODO: Verify if "all" strings should be enumerated - homomorphisms only?
  typedef std::function<void(const map_t &)> mappable_t;
  auto generateAllMaps = [&decomp, &nextMap](unsigned int p, mappable_t f)
  {
    auto length = decomp[p].size();
    std::vector<unsigned int> firstEnumeration(length, 0),
                              enumeration(firstEnumeration);

    auto bag = std::vector<unsigned int>(decomp[p].begin(), decomp[p].end());
    do
    {
      //Build homomorphism
      map_t homo;
      for(auto i = 0u; i < length; i++)
      {
        homo[bag.at(i)] = enumeration.at(i);
      }

      f(homo);
    } while(nextMap(enumeration));
  };

  for(auto p : ordering)
  {
    auto children = getChilds(p);
    //Start node (leafs)
    if(children.size() == 0)
    {
      assert(decomp[p].size() == 1);
      auto v = *decomp[p].begin();
      for(auto a = 0u; a < boost::num_vertices(host); ++a)
      {
        map_t key({std::make_pair(v,a)});
        counts[p][key] = 1;
      }
    }
    //Introduce node
    else if(children.size() == 1 &&
            decomp[p].size() == decomp[*children.begin()].size() + 1)
    {
      auto q = *children.begin();
      auto v = getUniqueSetDiff(p, q,children);

      //build Sq based on E(Gp)
      std::set<unsigned int> Sq;
      for(auto u : decomp[q])
      {
        //Check if there's any edge in the induced subgraph, involving u
        for(auto cand : decomp[p])
        {
          if(boost::edge(u, cand, pattern).second)
          {
            Sq.insert(u);
            break;
          }
        }
      }

      //forall f in Fq
      auto mappable = [p, q, v, &Sq, &host, &counts](map_t homo)
      {
        //for all a in V(H)
        for(auto a = 0u; a < boost::num_vertices(host); ++a)
        {
          //Determine which case applies
          bool allEdgesExist = true;
          for(auto u : Sq)
          {
            assert(homo.count(u) == 1);
            auto transformed = homo[u];
            if(!boost::edge(transformed, a, host).second)
            {
              allEdgesExist = false;
              break;
            }
          }

          auto extendedHomo = homo;
          extendedHomo[v] = a;
          if(allEdgesExist)
          {
            counts[p][extendedHomo] = counts[q][homo];
          }
          else
          {
            counts[p][extendedHomo] = 0;
          }
        }
      };
      generateAllMaps(p, mappable);

      counts.erase(q);
    }
    //Forget node
    else if(children.size() == 1 &&
            decomp[p].size() == decomp[*children.begin()].size() - 1)
    {
      auto q = *children.begin();
      auto v = getUniqueSetDiff(q, p, children);

      auto mappable = [v, p, q, &host, &counts](map_t homo)
      {
        //Sum over all host vertices
        unsigned int homoSum = 0;
        for(auto a = 0u; a < boost::num_vertices(host); ++a)
        {
          auto extendedHomo = homo;
          extendedHomo[v] = a;
          homoSum += counts[q][extendedHomo];
        }
        counts[p][homo] = homoSum;
      };
      generateAllMaps(p, mappable);
      counts.erase(q);
    }
    //Join node
    else
    {
      assert(children.size() == 2);
      auto it = children.begin();
      auto q1 = *it++;
      auto q2 = *it++;

      auto mappable = [p, q1, q2, &counts](const map_t & homo)
      {
        assert(counts[q1].count(homo) == 1 &&
               counts[q2].count(homo) == 1);

        counts[p][homo] = counts[q1][homo] * counts[q2][homo];
      };
      generateAllMaps(p, mappable);

      counts.erase(q1);
      counts.erase(q2);
    }

    visited.at(p) = true;
  }

  return counts[root][map_t()];
}

//TODO: Add implementation of counting injective homomorphisms

}
