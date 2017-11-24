# LibFastVC 
Fast iterative calculation of the minimal vertex cover.
Algorithm taken from [Shaowei Cai/FastVC](http://lcs.ios.ac.cn/~caisw/MVC.html) and converted to header only library.

The algorithm takes a graph in DIMACS format as input and calculates the
minimal vertex cover / independent set. During the calculation, approximations
are provided.

## Compiling

This version is header only and works with C++11. Just incluce `fastvc.hpp`.

## How to use

The standalone version uses LibFastVC and calculates the MVC of a graph given
as a text file in dimacs format.

```bash
make
# read file frb45-21-1.mis
# stop calculation only if time is up
# use 10 seconds for calculation
./fastvc frb45-21-1.mis 0 10
```

**Integration with BGL**

For integration with the Boost Graph Library (BGL) a converter is provided in
[this gist](https://gist.github.com/fmoessbauer/163b9928ae9170cfe2651173f416314b).

```cpp
// convert graph to dimacs
std::stringstream os;
write_dimacs(col_graph, os);

// initialize solver and solve
FastVC solver(os, boost::num_vertices(col_graph)-size_is, max_sec, verbose);
solver.cover_LS();
std::vector<int> solution = std::move(solver.get_independent_set());
```

## References

For standalone solvers see also http://lcs.ios.ac.cn/~caisw/MVC.html

