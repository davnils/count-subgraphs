#pragma once

#include <boost/numeric/ublas/matrix.hpp>
#include <map>
#include <set>
#include <vector>

#include "Utils.hpp"

namespace Count { namespace Basic {

typedef std::vector<std::map<std::set<unsigned int>, long>> parity_vec_t;
typedef std::pair<parity_vec_t, parity_vec_t> parity_t;

long calculateIntersectionTransform(
  const unsigned int,
  const std::set<unsigned int> &,
  const unsigned int,
  const unsigned int,
  std::function<long(const std::set<unsigned int> &)>
  );

long calculateParityTransform(
  const unsigned int,
  const std::set<unsigned int> &,
  const unsigned int,
  const unsigned int,
  const unsigned int,
  std::function<long(const std::set<unsigned int> &)>
  );

std::vector<std::vector<unsigned long long>> calculateAllTuples(
  const unsigned int,
  const unsigned int
  );

long getCoefficient(
  const unsigned int,
  const unsigned int,
  const unsigned int
  );

std::pair<parity_vec_t, parity_vec_t> buildTpValues(
  const Utils::partition_triple_t &,
  const unsigned int,
  const unsigned int,
  const double
  );

std::pair<bool, std::vector<long long>> buildSystem(
  const Utils::partition_triple_t &,
  const unsigned int,
  const unsigned int,
  const double
  );

/*matrix_t buildMatrixRepresentation(
  const std::vector<unsigned int> &,
  const std::vector<long> &,
  const unsigned int
  );*/

} }
