#include "TestTreeDecomposition.hpp"
#include "TreeDecomposition.hpp"
#include "TestHomomorphism.hpp"
#include "TestCount.hpp"
#include "TestBasicEquation.hpp"
#include "TestSimpleEquation.hpp"

int main()
{
  using namespace Count::Tree::Test;
  testBinaryDecomposition(std::cout);
  testNiceTreeDecomposition(std::cout);

  using namespace Count::Homomorphism::Test;
  testStingyOrdering(std::cout);
  testNaiveCounting(std::cout);
  testCountHomomorphisms(std::cout);
  testSimpleHomomorphism(std::cout);

  using namespace Count::Test;
  testCountInjective(std::cout);
  testPartitionGraph(std::cout);
  testCountInjectiveTriplesPath(std::cout);
  testCountInjectiveTriplesKPath(std::cout);
  testCountInjectiveTriples(std::cout);

  using namespace Count::Basic::Test;
  testLTable(std::cout);
  testBasicEquations(std::cout);

  using namespace Count::Simple::Test;
  testSymDiffProduct(std::cout);
  testSimpleEquations(std::cout);

  testCountIsoSubgraphs(std::cout);

  return 0;
}
