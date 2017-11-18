# Improved-NuMVC
An improved version, algorithm based on [khz1995/Improved-NuMVC](https://github.com/khz1995/Improved-NuMVC).

The algorithm takes a graph in DIMACS format as input and calculates the
minimal vertex cover / independent set. During the calculation, approximations
are provided.

## Compiling

This version is header only and works with C++11. Just incluce `tsewf.hpp`.

## How to use

**Integration with BGL**

For integration with the Boost Graph Library (BGL) a converter is provided in
`boost_graph.hpp`.

```cpp
// convert graph to dimacs
std::stringstream os;
write_dimacs(col_graph, os);

// initialize solver and solve
TSEWF solver(os, boost::num_vertices(col_graph)-num_bricks, max_sec, verbose);
solver.cover_LS();
std::vector<int> solution = std::move(solver.get_independent_set());
```

## References

For standalone solvers see also http://lcs.ios.ac.cn/~caisw/MVC.html

