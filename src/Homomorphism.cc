#include <boost/graph/adjacency_list.hpp>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "Homomorphism.hpp"
#include "InducedGraph.hpp"
#include "TreeDecomposition.hpp"
#include "Count.hpp"
#include "Utils.hpp"

namespace Count { namespace Homomorphism {

using namespace Utils;

/**
 * Calculate a stingy ordering satisfying critieria 1 and 2 (Diaz)
 * (1) the first element is a child, the last element is the root
 * (2) all parents occur after their childs
 *
 * only space bound affected by lack of critiera 3
 *
 * @param tree Tree to be used in calculations.
 * @param root Root node in the provided tree.
 * @return List of vertices, a stingy ordering.
 */
std::list<unsigned int> calculateStingyOrdering(
  const Tree::tree_decomp_t & tree,
  const unsigned int root
  )
{
  std::list<unsigned int> stingy;

  //BFS from root
  std::vector<bool> visited(boost::num_vertices(tree), false);
  std::queue<unsigned int> work;
  work.push(root);
  visited.at(root) = true;

  while(!work.empty())
  {
    auto v = work.front(); work.pop();
    stingy.push_front(v);

    for(auto it = boost::adjacent_vertices(v, tree); it.first != it.second; ++it.first)
    {
      auto u = *it.first;
      if(!visited.at(u))
      {
        visited.at(u) = true;
        work.push(u);
      }
    }
  }

  return stingy;
}

/**
 * Extract the subtree indicated by the root node and
 * some subroot node.
 *
 * @param tree Tree to be used.
 * @param root Root node in the input tree.
 * @param subRoot Subroot of tree to be extracted.
 */
std::set<unsigned int> extractSubTree(
  const Tree::tree_decomp_t & tree,
  const unsigned int root,
  const unsigned int subRoot
  )
{
  std::set<unsigned int> subTreeVertices;
  std::vector<bool> visited(boost::num_vertices(tree), false);
  std::queue<unsigned int> work;

  work.push(root);
  visited.at(root) = true;

  while(!work.empty())
  {
    auto v = work.front(); work.pop();

    if(v == subRoot)
    {
      work = decltype(work)();
      subTreeVertices.insert(v);
    }

    for(auto it = boost::adjacent_vertices(v, tree); it.first != it.second; ++it.first)
    {
      auto cand = *it.first;
      if(!visited.at(cand))
      {
        work.push(cand);
        visited.at(cand) = true;
      }
    }
  }

  return subTreeVertices;
}

/**
 * Count the number of homomorphisms between a pattern graph
 * and some host graph, extending the provided injective homomorphism.
 *
 * @param pattern Pattern graph (source).
 * @param decomp Decomposition of the provided pattern graph.
 * @param root Root node of the decomposition.
 * @param host Host graph (destination).
 * @param homoRef Provided injective homomorphism.
 * @return The number of valid extensions.
 */
unsigned long long countHomomorphisms(
  const Tree::undirected_graph_t & pattern,
  const Tree::tree_decomp_t & decomp,
  const unsigned int root,
  const Tree::undirected_graph_t & host,
  const map_t & homoRef
  )
{
  std::map<unsigned int, std::map<map_t, count_t>> counts;

  if(boost::num_vertices(host) == 0)
  {
    assert(boost::num_vertices(pattern) != 0);
    return 0;
  }

  auto ordering = calculateStingyOrdering(decomp, root);
  std::vector<bool> visited(ordering.size(), false);

  //TODO: Refactor
  auto getChilds = [&decomp, &visited](unsigned int v)
  {
    std::set<unsigned int> neighbours;
    for(auto it = boost::adjacent_vertices(v, decomp); it.first != it.second;
        ++it.first)
    {
      if(visited.at(*it.first))
      {
        neighbours.insert(*it.first);
      }
    }

    return neighbours;
  };

  auto nextMap =
    [&host]
    (std::vector<unsigned int> & vec)
  {
    auto const alphabetSize = boost::num_vertices(host);
    auto const length = vec.size();

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

  //TODO: Refactor
  auto getUniqueSetDiff =
    [&decomp]
    (unsigned int first, unsigned int second, const std::set<unsigned int> & children)
  {
    std::vector<unsigned int> diff(1);
    std::set_difference(decomp[first].begin(), decomp[first].end(),
                        decomp[second].begin(), decomp[second].end(),
                        diff.begin());
    return diff.at(0);
  };

  typedef std::function<void(const map_t &)> mappable_t;

  auto generateAllMaps =
    [&decomp, &nextMap, &homoRef]
    (const unsigned int p, mappable_t f)
  {
    auto const length = decomp[p].size();

    if(length == 0)
    {
      map_t homo;
      f(homo);
      return;
    }

    std::vector<unsigned int> enumeration(length, 0);

    auto isValidExtension =
      [&homoRef]
      (const map_t & cand)
    {
      for(auto const & entry : cand)
      {
        if(homoRef.count(entry.first) == 1 &&
           homoRef.at(entry.first) != entry.second)
        {
          return false;
        }
      }

      return true;
    };

    auto bag = std::vector<unsigned int>(decomp[p].begin(), decomp[p].end());
    do
    {
      //Construct map
      map_t homo;
      for(auto i = 0u; i < length; i++)
      {
        homo[bag.at(i)] = enumeration.at(i);
      }

      if(isValidExtension(homo))
      {
        f(homo);
      }

    } while(nextMap(enumeration));
  };

  for(auto const p : ordering)
  {
    auto children = getChilds(p);

    //Start node (leafs)
    if(children.size() == 0)
    {
      assert(decomp[p].size() == 1);
      auto v = *decomp[p].begin();

      //already included in homomorphism => one option
      auto homoIt = homoRef.find(v);
      if(homoIt != std::end(homoRef))
      {
        map_t key({*homoIt});
        counts[p][key] = 1;
      }
      //enumerate all host vertices otherwise
      else
      {
        for(auto a = 0u; a < boost::num_vertices(host); ++a)
        {
          map_t key({std::make_pair(v,a)});
          counts[p][key] = 1;
        }
      }
    }
    //Introduce node
    else if(children.size() == 1 &&
            decomp[p].size() == decomp[*children.begin()].size() + 1)
    {
      auto q = *children.begin();
      auto v = getUniqueSetDiff(p, q,children);

      auto subVertices = extractSubTree(decomp, root, p);
      assert(!subVertices.empty());

      //take the union of all bags from vertices in the subtree
      std::set<unsigned int> v_p;
      for(auto ext : subVertices)
      {
        v_p = Utils::mergeSets(v_p, decomp[ext]);
      }

      //S_q = all u s.t. (u,v) edge in G_p
      std::set<unsigned int> Sq;
      for(auto u : decomp[q])
      {
        //Check if there's any edge in the induced subgraph, involving u
        if(boost::edge(u, v, pattern).second &&
           v_p.count(u) == 1                 &&
           v_p.count(v) == 1)
        {
          Sq.insert(u);
        }
      }

      //forall f in Fq
      auto mappable = [p, q, v, &Sq, &host, &counts, &homoRef](map_t homo)
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
          
          //check if part of provided homomorphism
          auto allowed = true;
          if(homoRef.count(v) == 1 && homoRef.at(v) != a)
          {
            allowed = false;
          }

          auto extendedHomo = homo;
          extendedHomo[v] = a;
          if(allEdgesExist && allowed)
          {
            assert(counts[q].count(homo) == 1);
            counts[p][extendedHomo] = counts[q][homo];
          }
          else
          {
            counts[p][extendedHomo] = 0;
          }
        }
      };
      generateAllMaps(q, mappable);

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
          if(counts[q].count(extendedHomo))
          {
            homoSum += counts[q].at(extendedHomo);
          }
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

  return counts.at(root).at(map_t());
}

/**
 * Check if the provided homomorphism is valid.
 *
 * @param sourceGraph From this graph.
 * @param targetGraph To this graph.
 * @param homo A homomorphism.
 * @return True if valid.
 */
bool isValidPartialHomomorphism(
  const Tree::undirected_graph_t & sourceGraph,
  const Tree::undirected_graph_t & targetGraph,
  const map_t & homo
  )
{
  auto mappedVertices = extractKeys<unsigned int, unsigned int>(homo);

  for(auto const source : mappedVertices)
  {
    for(auto it = boost::adjacent_vertices(source, sourceGraph);
        it.first != it.second;
        ++it.first)
    {
      auto other = *it.first;
      if(mappedVertices.count(other) == 1 &&
         !boost::edge(homo.at(source), homo.at(other), targetGraph).second)
      {
        return false;
      }
    }
  }

  return true;
}

} }
