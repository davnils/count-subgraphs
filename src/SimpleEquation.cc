#include <algorithm>

#include "BasicEquation.hpp"
#include "SimpleEquation.hpp"

namespace Count { namespace Simple {

using namespace Utils;

/* TODO: REMOVE THESE */
std::ostream & operator<<(std::ostream & os, const std::vector<unsigned int> & vec)
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
std::map<std::set<unsigned int>, unsigned int> buildMatrixIndexMap(
  const unsigned int n,
  const unsigned int l
  )
{
  std::map<std::set<unsigned int>, unsigned int> output;

  assert(l <= n && l % 2 == 0);

  //position of the current set
  unsigned int indexCount = 0;

  auto evalEntry =
    [&indexCount, &output]
    (const std::vector<bool> & markers)
  {
    auto vec = extractSubset(markers);
    auto entry = std::set<unsigned int>(std::begin(vec), std::end(vec));
    output[entry] = indexCount++;
  };
  applyOnSubsets(n, l/2, evalEntry);

  return output;
}

/**
 *
 */
boost::numeric::ublas::matrix<unsigned long long> buildMatrix(
  const unsigned int n,
  const unsigned int rowSubsetParam,
  const unsigned int colSubsetParam,
  std::function<unsigned long long(std::set<unsigned int>)> eval
  )
{
  using namespace boost::numeric::ublas;

  auto totalRows = std::pow(n, rowSubsetParam);
  auto totalCols = std::pow(n, colSubsetParam);
  matrix<unsigned long long> output(totalRows, totalCols);
  output.assign(zero_matrix<unsigned long long>(totalRows, totalCols));

  unsigned int row = 0;
  auto evalRow =
    [&row, n, colSubsetParam, eval, &output]
    (const std::vector<bool> & markers)
  {
    auto rowVec = extractSubset(markers);
    auto rowSet = std::set<unsigned int>(std::begin(rowVec), std::end(rowVec));

    unsigned int col = 0;
    auto evalCol =
      [row, &col, eval, &output, &rowSet]
      (const std::vector<bool> & markers)
    {
      auto colVec = extractSubset(markers);
      auto colSet = std::set<unsigned int>(std::begin(colVec), std::end(colVec));

      if(buildSetIntersection<unsigned int>(rowSet, colSet).empty())
      {
        output(row, col) = eval(mergeSets(rowSet, colSet));
      }

      ++col;
    };
    applyOnSubsets(n, colSubsetParam, evalCol);

    ++row;
  };
  applyOnSubsets(n, rowSubsetParam, evalRow);

  return output;
}

/**
 *
 */
std::map<std::set<unsigned int>, long> calculateDiffTransform(
  const partition_triple_t & triple,
  const unsigned int q,
  const unsigned int n,
  const unsigned int minL,
  const unsigned int maxL
  )
{
  std::map<std::set<unsigned int>, long> result;

  //consider every l in the given range
  for(auto l = minL; l <= maxL; ++l)
  {
    if(l % 2 != 0)
    {
      continue;
    }

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

    //calculate matrix FG based on the current l-value
    auto matrixF = buildMatrix(n, l/2, q - l/2, fWrapper);
    auto matrixG = buildMatrix(n, q - l/2, l/2, gWrapper);
    auto matrixFG = prod(matrixF, matrixG);
    auto indexMap = buildMatrixIndexMap(n, l);

    //enumerate all subsets 'D'
    auto evalSubset =
      [l, q, n, &triple, &matrixFG, &indexMap, &result]
      (const std::vector<bool> & markers)
    {
      auto d = extractSubset(markers);
      auto dSet = std::set<unsigned int>(std::begin(d), std::end(d));

      auto iRange = d;
      long long sum = 0;

      //enumerate all subsets 'I'
      auto evalTerm =
        [&triple, &dSet, &matrixFG, &indexMap, &iRange, &sum]
        (const std::vector<bool> & markers)
      {
        auto i = extractSubset(markers, iRange);
        auto iSet = std::set<unsigned int>(std::begin(i), std::end(i));

        //read the corresponding entry in FG
        auto diffSet = buildSetDifference<unsigned int>(dSet, iSet);
        sum += matrixFG(indexMap.at(iSet), indexMap.at(diffSet));
      };
      applyOnSubsets(iRange.size(), l/2, evalTerm);

      result[dSet] = sum;
    };
    applyOnSubsets(n, l, evalSubset);
  }

  return result;
}

/**
 *
 */
std::vector<std::map<std::set<unsigned int>, long>> calculateIntersectionTransform(
  const partition_triple_t & triple,
  const unsigned int q,
  const unsigned int n,
  const unsigned int minL,
  const unsigned int maxL,
  const double gamma
  )
{
  auto s = std::floor((1 + 2*gamma)*q + 1);
  std::vector<std::map<std::set<unsigned int>, long>> result(s+1);

  for(auto l = minL; l <= maxL; ++l)
  {
    for(auto t = 0; t <= s; ++t)
    {
      auto evalSubset =
        [t, q, n, &triple, &result]
        (const std::vector<bool> & markers)
      {
        auto d = extractSubset(markers);
        auto dSet = std::set<unsigned int>(std::begin(d), std::end(d));

        auto hWrapper = [&triple](const std::set<unsigned int> & index)
        {
          if(triple.h.count(index) == 0)
          {
            return 0UL;
          }
          return triple.h.at(index);
        };

        result.at(t)[dSet] =
          Basic::calculateIntersectionTransform(n, dSet, q, t, hWrapper);
      };
      applyOnSubsets(n, l, evalSubset);
    }
  }

  return result;
}

/**
 *
 */
std::vector<long> solveIndeterminates(
  const partition_triple_t & triple,
  const unsigned int q,
  const unsigned int n,
  const unsigned int numSimple,
  const double gamma
  )
{
  //solve for all indeterminates
  std::vector<long> xVec;
  unsigned int lMin = n, lMax = 0;

  std::vector<unsigned int> indices;
  for(unsigned int j = 0u; j <= 3*q; ++j)
  {
    if(indices.size() == numSimple)
    {
      break;
    }

    if(j % 2 == q % 2)
    {
      indices.push_back(j);
      lMin = std::min(lMin, q - j);
      lMax = std::max(lMax, q + j);
    }
  }

  if(indices.empty())
  {
    return xVec;
  }

  //std::cerr << "lMin=" << lMin << ", lMax=" << lMax << ", q=" << q << ", n=" << n << std::endl;

  //std::cerr << "calculateDiffTransform()" << std::endl;
  auto diffTransform = calculateDiffTransform(triple, q, n, lMin, lMax);
  //std::cerr << "calculateIntersectionTransform()" << std::endl;
  auto intersectionTransform = calculateIntersectionTransform(triple, q, n, lMin, lMax, gamma);

  //solve each indeterminate
  for(auto j : indices)
  {
    //std::cerr << "j=" << j << std::endl;
    long total = 0;
    for(unsigned int l = q - j; l <= q + j; ++l)
    {
      //std::cerr << "l=" << l << std::endl;
      //enumerate all subsets D
      auto evalSubset =
        [q, l, j, &diffTransform, &intersectionTransform, &total]
        (const std::vector<bool> & markers)
      {
        auto d = extractSubset(markers);
        auto dSet = std::set<unsigned int>(std::begin(d), std::end(d));
        total += diffTransform.at(dSet) * intersectionTransform.at((q + l - j)/2).at(dSet);
      };
      //std::cerr << "generating subsets of size " << l << " over alphabet of size " << n << std::endl;
      applyOnSubsets(n, l, evalSubset);
    }

    xVec.push_back(total);
  }

  return xVec;
}


} }
