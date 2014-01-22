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

/**
 * Update the indices of an homomorphism based on the provided maps.
 *
 * @param prev The old homomorphism.
 * @param indices Map from old to new vertex, or -1 if excluded.
 * @param sourceMap Flag indicating if the update is performed on
 *                  source (true) or target (false) vertices.
 * @return Homomorphism with updated indices.
 */
static inj_homo_t updateHomoIndexing(
  const inj_homo_t & prev,
  const std::vector<int> & indices,
  bool sourceMap
  )
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


/**
 * Count the number of injective homomorphisms extending the supplied inj. homo.,
 * for all subsets of size |pattern| - |assigned|.
 * Corresponds to the result by Fomin.
 *
 * @param pattern Pattern graph (source).
 * @param host Host graph (target).
 * @param excluded Set of excluded vertices (drawn from host graph).
 * @param homo Existing injective homomorphism.
 * @return Map from subsets to the corresponding count.
 */
std::map<std::set<unsigned int>, unsigned long> countInjective(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host,
  const std::set<unsigned int> & excluded,
  const inj_homo_t & homo
  )
{
  //build tree decomposition of pattern graph
  auto decomp = Tree::convertToNiceDecomposition(Tree::buildTreeDecomposition(pattern));

  std::set<unsigned int> homoMappedSubset;
  for(auto const & entry : homo)
  {
    homoMappedSubset.insert(entry.second);
  }

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

    auto result = Homomorphism::countHomomorphisms(pattern,
                                                   decomp.first,
                                                   decomp.second,
                                                   subResult.first,
                                                   subHomo);

    homoCount[vertexSet] = result;
  };

  //size of the largest subset (later combined with existing homomorphism)
  auto const maxLength = boost::num_vertices(pattern) - homoMappedSubset.size();

  //number of available vertices in the host graph
  auto const hostSetSize = freeVertices.size();

  for(auto length = 0u; length <= maxLength; ++length)
  {
    applyOnSubsets(hostSetSize, length, evaluateHomomorphisms);
  }

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
      long coeff = std::pow(-1, (long)(maxLength - subsetVertices.size()));
      auto term = homoCount[std::set<unsigned int>(std::begin(subsetVertices),
                                                   std::end(subsetVertices))];
      count += coeff * (long)term;
    };

    //enumarate all vertex subsets of 'set'
    for(auto length = 0u; length <= vertices.size(); ++length)
    {
      applyOnSubsets(vertices.size(), length, f);
    }

    auto vertexSet = std::set<unsigned int>(std::begin(vertices), std::end(vertices));

    //take the provided injective homomorphism into account
    if(vertexSet.empty())
    {
      count = 1;
    }

    assert(count >= 0);
    injCount[vertexSet] = (unsigned long)count;
  };

  //enumerate all subsets Q of V(H) \ phi(assigned), of size exactly |P| - |a|
  applyOnSubsets(hostSetSize, maxLength, evaluateZetaTransform);

  //return the total number of injective homomorphisms
  return injCount;
}

/**
 * Count the number of injective homomorphisms between different subgraphs.
 *
 * @param pattern Pattern graph (source).
 * @param host Host graph (target).
 * @return Struct containing counts and partitioning.
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
 * Count the number of injective homomorphisms between different subgraphs.
 *
 * @param pattern Pattern graph (source).
 * @param host Host graph (target).
 * @param givenPartition Partition to be used.
 * @return Struct containing counts and partitioning.
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
      //only keep vertices from S
      auto patternFilter = mergeSets(subResults.part.L, subResults.part.S);

      //consider the pattern subgraph
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      subResults.counts[homo].f =
        countInjective(subPattern.first,
                       host,
                       mapVertices(homo, subResults.part.T),
                       filterHomo(subHomo, updateIndices(subPattern.second, subResults.part.S)));
    }

    //P[S `union` M `union` T] -> H[phi(S) `union` B `union` phi(T)]
    {
      //consider the pattern subgraph
      auto patternFilter = mergeSets(subResults.part.S, mergeSets(subResults.part.M, subResults.part.T));
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      std::set<unsigned int> emptySet;
      subResults.counts[homo].g = countInjective(subPattern.first, host, emptySet, subHomo);
    }

    //P[R `union` T] -> H[C `union` phi(T)]
    {
      //only keep vertices from T
      auto patternFilter = mergeSets(subResults.part.R, subResults.part.T);

      //consider the pattern subgraph
      auto subPattern = InducedGraph::buildInducedSubGraphTagged<Tree::undirected_graph_t>(pattern, patternFilter);
      auto subHomo = updateHomoIndexing(homo, subPattern.second, true);
      subResults.counts[homo].h =
        countInjective(subPattern.first,
                       host,
                       mapVertices(homo, subResults.part.S),
                       filterHomo(subHomo, updateIndices(subPattern.second, subResults.part.T)));
    }
  };

  //generate all injective homormorphisms mapping S and T
  //take all |S|+|T|-sized subsets of V(H) and enumerate all vertex permutations for each subset
  auto evaluateCandidate =
    [&evaluateInjectiveHomomorphism, &host, &pattern, &subResults]
    (const std::vector<bool> & candidateIndices)
  {
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
 * Calculate an entry in the A matrix based on row and column.
 *
 * @param n Total dimension.
 * @param i Row.
 * @param j Column.
 * @return The corresponding coefficient.
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
 * Solve an integer linear equation system on the format Ax=y.
 *
 * @param A A matrix.
 * @param y y vector.
 * @return Vector containing solutions to all indeterminates.
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

/**
 * Build a linear equation system based on triples calculated beforehand.
 *
 * @param n Number of host vertices.
 * @param q Size of subsets being surveyed, typically |P| = 3q.
 * @param triple Triples containing injective counts for various subgraphs.
 * @return A and y matrices, from Ax = y.
 */
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

  //build matrices from basic equations
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

  //optionally, solve indeterminates and add them to the system
  auto solved = Simple::solveIndeterminates(triple, q, n, e - d, gamma);
  for(auto i = 0u; i < e - d; ++i)
  {
    A(d + i, i) = 1;
    y(d + i, 0) = solved.at(i);
  }

  return std::make_pair(A, y);
}

/**
 * Count the number of subgraphs isomorphic to the pattern graph.
 * Limited to pattern graphs admitting an even 5-partitioning, 
 * where the size of L, M, and R, equals.
 *
 * @param pattern bla.
 * @param host bla.
 * @return Boolean flag indicating success and the count.
 */
std::pair<bool, unsigned long long> countIsoSubgraphs(
  const Tree::undirected_graph_t & pattern,
  const Tree::undirected_graph_t & host
  )
{
  //parameter controlling the balance between basic and simple equations
  const double gamma = 0.5;

  //build 5-partitioning
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
