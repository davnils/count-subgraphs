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
 * calculate a stingy ordering satisfying critieria 1 and 2 (Diaz)
 * (1) the first element is a child, the last element is the root
 * (2) all parents occur after their childs
 *
 * only space bound affected by lack of critiera 3
 */
std::list<unsigned int> calculateStingyOrdering(
  const Tree::tree_decomp_t & tree,
  const unsigned int root)
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
 *
 */
std::set<unsigned int> extractSubTree(
  const Tree::tree_decomp_t & tree,
  const unsigned int root,
  const unsigned int subRoot)
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

//TODO: remove these
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
unsigned long long countHomomorphisms(
  const Tree::undirected_graph_t & pattern,
  const Tree::tree_decomp_t & decomp,
  const unsigned int root,
  const Tree::undirected_graph_t & host,
  const map_t & homoRef)
{
  std::map<unsigned int, std::map<map_t, count_t>> counts;

  if(boost::num_vertices(host) == 0)
  {
    assert(boost::num_vertices(pattern) != 0);
    //std::cerr << "ignoring empty host graph" << std::endl;
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

    //std::cerr << "alphabetSize=" << alphabetSize  << ", length=" << length << std::endl;

    auto nextIndex = 0;
    while(nextIndex < length &&
          vec.at(nextIndex) == alphabetSize - 1)
    {
      vec.at(nextIndex) = 0;
      ++nextIndex;
    }

    if(nextIndex == length)
    {
      //std::cerr << "returning false" << std::endl;
      return false;
    }

    vec.at(nextIndex)++;

    //std::cerr << "returning true" << std::endl;
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
    //std::cerr << "generateAllMaps(p=" << p << ")" << std::endl;
    auto const length = decomp[p].size();
    //std::cerr << "length=" << length << std::endl;

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

      //TODO: Updated from continue
      if(isValidExtension(homo))
      {
        f(homo);
      }

    } while(nextMap(enumeration));
  };


  //TODO: What about empty host graphs?

  for(auto const p : ordering)
  {
    //std::cerr << "p=" << p << std::endl;
    auto children = getChilds(p);

    //Start node (leafs)
    if(children.size() == 0)
    {
      //std::cerr << ">>>> START" << std::endl;
      assert(decomp[p].size() == 1);
      auto v = *decomp[p].begin();

      //already included in homomorphism => one option
      auto homoIt = homoRef.find(v);
      if(homoIt != std::end(homoRef))
      {
        //std::cerr << "v=" << v << " already in homo." << std::endl;
        map_t key({*homoIt});
        counts[p][key] = 1;
      }
      //enumerate all host vertices otherwise
      else
      {
        //std::cerr << "enumerating all host vertices with v=" << v << std::endl;
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
      //std::cerr << ">>>> INTRODUCE" << std::endl;
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

      //std::cerr << "S_q = " << Sq << std::endl;

      //forall f in Fq
      auto mappable = [p, q, v, &Sq, &host, &counts, &homoRef](map_t homo)
      {
        //std::cerr << "new homomorphism: " << homo << std::endl;
        //for all a in V(H)
        for(auto a = 0u; a < boost::num_vertices(host); ++a)
        {
          //std::cerr << "considering mapping " << v << " onto " << a << std::endl;
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
          
          //check if part of provided homomorphism (TODO: rewrite)
          auto allowed = true;
          if(homoRef.count(v) == 1 && homoRef.at(v) != a)
          {
            //std::cerr << "rejected" << std::endl;
            allowed = false;
          }

          auto extendedHomo = homo;
          extendedHomo[v] = a;
          if(allEdgesExist && allowed)
          {
            /*std::cerr << "writing counts[" << p << "] = counts[" << q << "][" << homo << "]"
                      << "(" << counts[q][homo] << ")" << std::endl;*/
            assert(counts[q].count(homo) == 1);
            counts[p][extendedHomo] = counts[q][homo];
          }
          else
          {
            //std::cerr << "writing 0" << std::endl;
            counts[p][extendedHomo] = 0;
          }
        }
      };
      generateAllMaps(q, mappable);

      counts.erase(q);
      //std::cerr << "end of intro" << std::endl;
    }
    //Forget node
    else if(children.size() == 1 &&
            decomp[p].size() == decomp[*children.begin()].size() - 1)
    {
      //std::cerr << ">>>> FORGET" << std::endl;
      auto q = *children.begin();
      auto v = getUniqueSetDiff(q, p, children);
      //std::cerr << "v=" << v << std::endl;

      auto mappable = [v, p, q, &host, &counts](map_t homo)
      {
        //std::cerr << "new homomorphism: " << homo << std::endl;
        //Sum over all host vertices
        unsigned int homoSum = 0;
        for(auto a = 0u; a < boost::num_vertices(host); ++a)
        {
          auto extendedHomo = homo;
          extendedHomo[v] = a;
          //std::cerr << "extending homo to: " << extendedHomo << std::endl;
          if(counts[q].count(extendedHomo))
          {
            //std::cerr << "adding " << counts[q].at(extendedHomo) << std::endl;
            homoSum += counts[q].at(extendedHomo);
          }
        }

        //std::cerr << "writing counts[" << p << "][" << homo << "]=" << homoSum << std::endl;
        counts[p][homo] = homoSum;
      };

      generateAllMaps(p, mappable);
      counts.erase(q);
      //std::cerr << "end of forget" << std::endl;
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

bool isValidPartialHomomorphism(
  const Tree::undirected_graph_t & sourceGraph,
  const Tree::undirected_graph_t & targetGraph,
  const map_t & homo)
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
