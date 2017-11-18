#ifndef BGL_WRITE_DIMACS_INC
#define BGL_WRITE_DIMACS_INC

/**
 * export boost graph to dimacs graph format
 */
template<typename graph_t>
void write_dimacs(const graph_t & graph, std::ostream & out){
  auto num_edges    = boost::num_edges(graph);
  auto num_vertices = boost::num_vertices(graph);
  boost::graph_traits<graph_t>::edge_iterator ei, ei_end;

  out << "p edge " << num_vertices << " " << num_edges << std::endl;
  for (tie(ei, ei_end) = boost::edges(graph); ei != ei_end; ++ei){
    auto source = boost::source(*ei, graph);
    auto target = boost::target(*ei, graph);
    auto sidx   = boost::get(boost::vertex_index, graph, source);
    auto tidx   = boost::get(boost::vertex_index, graph, target);
    out << "e " << sidx+1 << " " << tidx+1 << std::endl;
  }
}
#endif

