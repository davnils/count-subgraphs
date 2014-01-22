#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <map>
#include <set>

#include "BasicEquation.hpp"
#include "Count.hpp"
#include "Homomorphism.hpp"
#include "InducedGraph.hpp"
#include "Naive.hpp"
#include "SimpleEquation.hpp"
#include "TreeDecomposition.hpp"
#include "Utils.hpp"

namespace Count {

using namespace Utils;

std::ostream & operator<<(std::ostream & os, const std::vector<bool> & vec)
{
  os << " {";
  for(auto e : vec)
  {
    os << (char)('0' + e) << " ";
  }
  os << "} ";

  return os;
};

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
static inj_homo_t updateHomoIndexing(
  const inj_homo_t & prev, const std::vector<int> & indices, bool sourceMap)
{
  //update homomorphism with the new subgraph indexing, skip any excluded vertices
  inj_homo_t updated; 

  for(auto entry : prev)
  {
    int target;
    if(sourceMap)
    {
      if((target = indices.at(entry.first)) != -1)
      {
        updated[target] = entry.second;
      }
    }
    else
    {
      if((target = indices.at(entry.second)) != -1)
      {
        updated[entry.first] = target;
      }
    }
  }

  return updated;
}

//bool debug = false;


/**
 * counts the number of injective homomorphisms extending the supplied inj. homo.,
 * for all subsets of size |pattern| - |assigned|
 */
std::map<std::set<unsigned int>, unsigned long> countInjective(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host,
  const std::set<unsigned int> & excluded,
  const inj_homo_t & homo)
{
  /*std::cerr << "considering a given map of a candidate phi" << std::endl;
  std::cerr << "|P| = " << boost::num_vertices(pattern) << ", |H| = " << boost::num_vertices(host)
            << ", provided homomorphism: " << homo << std::endl;*/

  //Build tree decomposition of pattern graph
  auto decomp = Tree::convertToNiceDecomposition(Tree::buildTreeDecomposition(pattern));

  /*if(debug)
  {
    std::cerr << "------------------------------------------" << std::endl;
    std::cerr << "pattern graph" << std::endl;
    visualizeGraph(std::cerr, pattern);
    std::cerr << "host graph" << std::endl;
    visualizeGraph(std::cerr, host);
    std::cerr << "generated decomposition" << std::endl;
    Tree::visualizeDecomposition(std::cerr, decomp.first);
    std::cerr << "homo=" << homo << std::endl;
  }*/

  std::set<unsigned int> homoMappedSubset;
  for(auto const & entry : homo)
  {
    homoMappedSubset.insert(entry.second);
  }

  //std::cerr << "size of homoMapped: " << homoMappedSubset.size() << std::endl;

  //extract all host vertices not yet included
  std::vector<unsigned int> freeVertices;
  for(unsigned int cand = 0; cand < boost::num_vertices(host); ++cand)
  {
    if(homoMappedSubset.count(cand) == 0 && excluded.count(cand) == 0)
    {
      freeVertices.push_back(cand);
    }
  }

  //enumerate all subsets Q of V(H) \ phi(assigned), of size at most |P| - |a|
  //evaluate hom_phi(P, G[Q `union` phi(assigned)]) over all suitable Q
  std::map<std::set<unsigned int>, unsigned long> homoCount;
  auto evaluateHomomorphisms =
    [&pattern, &decomp, &host, &homo, &homoCount, &freeVertices, &homoMappedSubset]
    (const std::vector<bool> & subset)
  {
    auto subGraphVertices = extractSubset(subset, freeVertices);
    auto vertexSet = std::set<unsigned int>(std::begin(subGraphVertices), std::end(subGraphVertices));

    auto subResult = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>
                       (host, mergeSets(vertexSet, homoMappedSubset));

    auto subHomo = updateHomoIndexing(homo, subResult.second, false);

    //std::cerr << ">>>>>>>>>>>>>>>>>>>> considering a single homomorphism " << std::endl;
    //std::cerr << "homo: " << homo << std::endl;
    //std::cerr << "subHomo: " << subHomo << std::endl;

    /*std::cerr << "HOST GRAPH" << std::endl;
    visualizeGraph(std::cerr, subResult.first);
    std::cerr << "PATTERN GRAPH" << std::endl;
    visualizeGraph(std::cerr, pattern);
    std::cerr << "DECOMPOSITION GRAPH (r=" << decomp.second << ")" << std::endl;
    Tree::visualizeDecomposition(std::cerr, decomp.first);*/

    //std::cerr << "calling count()" << std::endl;
    auto result = Homomorphism::countHomomorphisms(pattern,
                                                   decomp.first,
                                                   decomp.second,
                                                   subResult.first,
                                                   subHomo);

    /*std::cerr << "(homo) writing injcount[" << vertexSet << "] = " << result
              << std::endl;*/
    homoCount[vertexSet] = result;
  };

  //size of the largest subset (later combined with existing homomorphism)
  auto const maxLength = boost::num_vertices(pattern) - homoMappedSubset.size();

  //number of available vertices in the host graph
  auto const hostSetSize = freeVertices.size();

  //std::cerr << ">>>>> iterating over all subsets of size at most " << maxLength << std::endl;
  for(auto length = 0u; length <= maxLength; ++length)
  {
    //std::cerr << "evaluating length=" << length << std::endl;
    applyOnSubsets(hostSetSize, length, evaluateHomomorphisms);
  }

  //std::cerr << "entering the second phase" << std::endl;

  //accumulate all values (using naive zeta transform), with -1 coefficient applied
  std::map<std::set<unsigned int>, unsigned long> injCount;
  auto evaluateZetaTransform =
    [&injCount, &homoCount, &freeVertices, &maxLength]
    (const std::vector<bool> & set)
  {
    //extract vertices associated with 'set'
    auto vertices = extractSubset(set, freeVertices);
    long long count = 0;

    //extract vertices corresponding to the subset and accumulate the corresponding
    //number of homomorphisms
    auto f =
      [&count, &vertices, &homoCount, &maxLength]
      (const std::vector<bool> & subset)
    {
      auto subsetVertices = extractSubset(subset, vertices);
      //std::cerr << "coeff= -1 ^ " << maxLength << " - " << subsetVertices.size() << std::endl; 
      long coeff = std::pow(-1, (long)(maxLength - subsetVertices.size()));
      auto term = homoCount[std::set<unsigned int>(std::begin(subsetVertices),
                                                   std::end(subsetVertices))];
      /*std::cerr << "adding subset: "
                << std::set<unsigned int>(std::begin(subsetVertices), std::end(subsetVertices))
                << "(" << coeff << " * " << term << ")" << std::endl;*/
      count += coeff * (long)term;
    };

    //enumarate all vertex subsets of 'set'
    for(auto length = 0u; length <= vertices.size(); ++length)
    {
      applyOnSubsets(vertices.size(), length, f);
    }

    auto vertexSet = std::set<unsigned int>(std::begin(vertices), std::end(vertices));

    //take the provided injective homomorphism into account <- should it be here? TODO
    if(vertexSet.empty())
    {
      count = 1;
    }

    //std::cerr << "injcount[" << vertexSet << "] = " << count << std::endl;
    assert(count >= 0);
    injCount[vertexSet] = (unsigned long)count;
  };

  //std::cerr << ">>>>>>> enumerating all full-sized subsets" << std::endl;
  //enumerate all subsets Q of V(H) \ phi(assigned), of size exactly |P| - |a|
  applyOnSubsets(hostSetSize, maxLength, evaluateZetaTransform);

  //return the total number of injective homomorphisms
  return injCount;
}

/**
 *
 */
subgraph_result countSubgraphTriples(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host
  )
{

  PatternPartition partition(pattern);
  return countSubgraphTriples(pattern, host, partition);
}

/**
 *
 */
subgraph_result countSubgraphTriples(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host,
  const PatternPartition & givenPartition
  )
{
  subgraph_result subResults(givenPartition);

  if(!subResults.part.partitionExists())
  {
    std::cerr << "No partition found, skipping." << std::endl;
    return subResults;
  }

  //evaluate a valid injective homo. mapping S,T -> Host
  auto evaluateInjectiveHomomorphism =
    [&subResults, &pattern, &host]
    (const inj_homo_t & homo)
  {
    auto filterHomo =
      []
      (const inj_homo_t & homo, const std::set<unsigned int> & include)
    {
      auto homoCopy = homo;
      auto it = std::begin(homoCopy);
      while(it != std::end(homoCopy))
      {
        if(include.count(it->first) == 0)
        {
          it = homoCopy.erase(it);
        }
        else
        {
          ++it;
        }
      }

      return homoCopy;
    };

    auto mapVertices =
      []
      (const inj_homo_t & homo, const std::set<unsigned int> & vertices)
    {
      std::set<unsigned int> mapped;

      for(auto v : vertices)
      {
        mapped.insert(homo.at(v));
      }

      return mapped;
    };

    auto updateIndices =
      []
      (const std::vector<int> & map, const std::set<unsigned int> & vertices)
    {
      std::set<unsigned int> mapped;

      for(auto v : vertices)
      {
        assert(map.at(v) != -1);
        mapped.insert(map.at(v));
      }

      return mapped;
    };

    //P[L `union` S] -> H[A `union` phi(S)]
    {
      /*std::cerr << "P[L `union` S] -> H[A `union` phi(S)] with s="
                << homo.at(*partition.S.begin()) << " and t="
                << homo.at(*partition.T.begin()) << std::endl;*/
      //Only keep vertices from S
      auto patternFilter = mergeSets(subResults.part.L, subResults.part.S);

      //Consider the pattern subgraph
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      subResults.counts[homo].f = countInjective(subPattern.first, host, mapVertices(homo, subResults.part.T), filterHomo(subHomo, updateIndices(subPattern.second, subResults.part.S)));
    }

    //P[S `union` M `union` T] -> H[phi(S) `union` B `union` phi(T)]
    {
      /*std::cerr << "P[S `union` M `union` T] -> H[phi(S) `union` B `union` phi(T)] with s="
                << homo.at(*partition.S.begin()) << " (" << partition.S << ") and m="
                << mVar << " and t="
                << homo.at(*partition.T.begin()) << " (" << partition.T << ")" << std::endl;*/
      //Consider the pattern subgraph
      auto patternFilter = mergeSets(subResults.part.S, mergeSets(subResults.part.M, subResults.part.T));
      //std::cerr << partition.S << " " << partition.M << " " << partition.T << std::endl;
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      std::set<unsigned int> emptySet;

      /*if(homo.count('c' - 'a') && homo.count('f' - 'a') && homo.count('b' - 'a') &&
         homo.at('c' - 'a') == ('d' - 'a') && homo.at('f' - 'a') == ('f' - 'a') && homo.at('b' - 'a') == ('b' - 'a'))
      {
        debug = true;
      }*/

      subResults.counts[homo].g = countInjective(subPattern.first, host, emptySet, subHomo);

      //debug = false;
    }

    //P[R `union` T] -> H[C `union` phi(T)]
    {
      /*std::cerr << "P[R `union` T] -> H[C `union` phi(T)] with phi(s)="
                << homo.at(*partition.S.begin()) << " and r="
                << rVar << " and phi(t)="
                << homo.at(*partition.T.begin()) << std::endl;*/
      //Only keep vertices from T
      auto patternFilter = mergeSets(subResults.part.R, subResults.part.T);

      //Consider the pattern subgraph
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      subResults.counts[homo].h = countInjective(subPattern.first, host, mapVertices(homo, subResults.part.S), filterHomo(subHomo, updateIndices(subPattern.second, subResults.part.T)));
    }
  };

  //generate all injective homormorphisms mapping S and T
  //take all |S|+|T|-sized subsets of V(H) and enumerate all vertex permutations for each subset
  auto evaluateCandidate =
    [&evaluateInjectiveHomomorphism, &host, &pattern, &subResults]
    (const std::vector<bool> & candidateIndices)
  {
    //std::cerr << "evaluateCandidate()" << std::endl;
    auto hostVertices = extractSubset(candidateIndices);

    std::vector<unsigned int> patternVertices(std::begin(subResults.part.S), std::end(subResults.part.S));
    patternVertices.insert(std::end(patternVertices), std::begin(subResults.part.T), std::end(subResults.part.T));

    //consider all permutations of this subset
    std::sort(std::begin(hostVertices), std::end(hostVertices));
    
    do
    {
      //extract the homomorphism
      inj_homo_t homo;
      for(auto i = 0u; i < hostVertices.size(); ++i)
      {
        homo[patternVertices.at(i)] = hostVertices.at(i);
      }

      //check if it is a valid (injective) homomorphism
      //the injective property holds by construction
      if(!Homomorphism::isValidPartialHomomorphism(pattern, host, homo))
      {
        continue;
      }

      //count
      evaluateInjectiveHomomorphism(homo);
    } while(std::next_permutation(std::begin(hostVertices), std::end(hostVertices)));
  };

  //enumerate all possible injective homomorphisms
  applyOnSubsets(boost::num_vertices(host),
                 subResults.part.S.size() + subResults.part.T.size(),
                 evaluateCandidate);

  return subResults;
}

/**
 *
 */
long getCoefficient(
  const unsigned int n,
  const unsigned int i,
  const unsigned int j
  )
{
  return std::pow((long)n - (long)j*2, (long)i);
}

/**
 *
 */
std::vector<int> solveLinearSystem(
  const boost::numeric::ublas::matrix<int> & A,
  const boost::numeric::ublas::matrix<int> & y
  )
{
  using namespace boost::numeric::ublas;

  matrix<double> Acopy(A);

  matrix<double> inverse(A.size1(), A.size1());
  inverse.assign(identity_matrix<double>(A.size1()));

  permutation_matrix<size_t> perm(A.size1());

  assert(lu_factorize(Acopy, perm) == 0);
  lu_substitute(Acopy, perm, inverse);

  auto tmp = prod(inverse, y);

  std::vector<int> solutions(A.size2());
  for(auto j = 0u; j < A.size2(); ++j)
  {
    solutions.at(j) = std::round(tmp(j, 0));
  }

  return solutions;
}

void printMatrix(std::ostream & os, const matrix_t & m)
{
  for(auto i = 0u; i < m.size1(); ++i)
  {
    for(auto j = 0u; j < m.size2(); ++j)
    {
      os << m(i, j) << " ";
    }
    os << std::endl;
  }
}

std::pair<matrix_t, matrix_t> buildSystem(
  const unsigned int n,
  const unsigned int q,
  const partition_triple_t & triple,
  const double gamma
  )
{
  auto equations = Basic::buildSystem(triple, q, n, gamma);
  auto const e = std::floor(3.0f*q/2.0f) + 1;
  auto const d = std::floor((3.0f/2.0f - gamma)*q) + 1;
  boost::numeric::ublas::matrix<int> A(e, e),
                                     y(e, 1);

  A.assign(boost::numeric::ublas::zero_matrix<int>(e, e));

  //build vector containing all indeterminates indices
  std::vector<unsigned int> indeterminates;
  for(auto j = 0u; indeterminates.size() < e; ++j)
  {
    if(j % 2 == q % 2)
    {
      indeterminates.push_back(j);
    }
  }

  for(auto i = 0u; i < d; ++i)
  {
    for(auto j = 0u; j < indeterminates.size(); ++j)
    {
      auto var = indeterminates.at(j);
      auto coeff = Count::getCoefficient(n, i, var);
      A(i, j) = coeff;
    }
    y(i, 0) = equations.second.at(i);
  }

  auto solved = Simple::solveIndeterminates(triple, q, n, e - d, gamma);
  for(auto i = 0u; i < e - d; ++i)
  {
    A(d + i, i) = 1;
    y(d + i, 0) = solved.at(i);
  }

  /*std::cerr << ">>> A:" << std::endl;
  printMatrix(std::cerr, A);
  std::cerr << ">>> y:" << std::endl;
  printMatrix(std::cerr, y);*/

  return std::make_pair(A, y);
}

std::pair<bool, unsigned long long> countIsoSubgraphs(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host
  )
{
  //parameter controlling the balance between basic and simple equations
  const double gamma = 0.5;

  //build 5-partition
  auto predicate = 
    []
    (PatternPartition * p)
  {
      return p->S.size() > 0 &&  p->T.size() > 0 &&
             p->L.size() == p->M.size() && p->L.size() == p->R.size();
  };

  PatternPartition part(pattern, predicate);

  if(!part.partitionExists())
  {
    return std::make_pair(false, 0);
  }

  auto n = boost::num_vertices(host);
  auto q = part.L.size();

  //construct all injective homomorphisms involving S and T
  auto injResults = countSubgraphTriples(pattern, host, part);

  unsigned long long injective = 0;
  for(auto const & triple : injResults.counts)
  {
    //add up the contribution gathered by solving the system
    auto system = buildSystem(n, q, triple.second, gamma);
    auto contrib = Count::solveLinearSystem(system.first, system.second).back();
    injective += contrib;
  }

  //calculate the number of automorphisms
  auto autoResults = countSubgraphTriples(pattern, pattern, part);
  auto automorphisms = Naive::evaluateDisjointTriples(autoResults.counts, pattern, part);

  assert(injective % automorphisms == 0);
  return std::make_pair(true, injective / automorphisms);
}

}
