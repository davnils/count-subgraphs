#pragma once

#include <map>
#include <set>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>

#include "Utils.hpp"

namespace Count { namespace Simple {

std::map<std::set<unsigned int>, unsigned int> buildMatrixIndexMap(
  const unsigned int,
  const unsigned int
  );

boost::numeric::ublas::matrix<unsigned long long> buildMatrix(
  const unsigned int,
  const unsigned int,
  const unsigned int,
  std::function<unsigned long long(std::set<unsigned int>)>
  );

std::map<std::set<unsigned int>, long> calculateDiffTransform(
  const Utils::partition_triple_t &,
  const unsigned int,
  const unsigned int,
  const unsigned int,
  const unsigned int
  );

std::vector<std::map<std::set<unsigned int>, long>> calculateIntersectionTransform(
  const Utils::partition_triple_t &,
  const unsigned int,
  const unsigned int,
  const unsigned int,
  const unsigned int,
  const double
  );

std::vector<long> solveIndeterminates(
  const Utils::partition_triple_t &,
  const unsigned int,
  const unsigned int,
  const unsigned int,
  const double
  );

} }
