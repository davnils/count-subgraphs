#include <boost/graph/graphviz.hpp>
#include <dlib/graph_utils.h>
#include <dlib/graph.h>

#include "TreeDecomposition.hpp"

namespace count {

class VertexSetWriter
{
public:
  VertexSetWriter(const tree_decomp_t & t)
    : m_tree(t)
  {
  }

  template <typename Item>
  void operator()(std::ostream & os, const Item & v) const
  {
    os << "[label=\"{";
    bool first = true;
    for(auto u : m_tree[v])
    {
      if(!first)
      {
        os << ",";
      }
      first = false;

      os << (char)('a' + u);
    }
    os <<  "}\"]";
  }

private:
  const tree_decomp_t & m_tree;
};

tree_decomp_t buildTreeDecomposition(const undirected_graph_t & inputGraph)
{
  dlib::graph<int>::kernel_1a_c inputTransformed;
  inputTransformed.set_number_of_nodes(boost::num_vertices(inputGraph));

  for(auto it = boost::edges(inputGraph); it.first != it.second; it.first++)
  {
    inputTransformed.add_edge(boost::source(*it.first, inputGraph),
                              boost::target(*it.first, inputGraph));
  }

  typedef dlib::set<unsigned long>::compare_1b_c set_t;
  dlib::graph<set_t, set_t>::kernel_1a_c internalDecomp;

  dlib::create_join_tree(inputTransformed, internalDecomp);

  tree_decomp_t outputDecomp(internalDecomp.number_of_nodes());
  for(auto v = 0; v < internalDecomp.number_of_nodes(); v++)
  {
    auto & it = internalDecomp.node(v).data;
    it.reset();
    while(it.move_next())
    {
      outputDecomp[v].insert(it.element());
    }

    for(auto u = v; u < internalDecomp.number_of_nodes(); u++)
    {
      if(internalDecomp.has_edge(v, u))
      {
        boost::add_edge(v, u, outputDecomp);
      }
    }
  }

  return outputDecomp;
}

void visualizeDecomposition(std::ostream & os, const tree_decomp_t & tree)
{
  boost::write_graphviz(os, tree, VertexSetWriter(tree));
}

}
