# LibMVC

[![CircleCI](https://circleci.com/gh/fmoessbauer/LibMVC.svg?style=shield)](https://circleci.com/gh/fmoessbauer/LibMVC)

LibMVC is a collection of fast iterative minimum vertex cover solvers.
Currently the following algorithms are implemented:

- NuMVC
- FastVC

The solvers take a graph in DIMACS format as input and calculate the
minimum vertex cover / independent set. During the calculation, approximations
are provided.

### Parallel Solver Adapter

The library includes an adapter for parallelizing any LibMVC interface
compatible solver:

```cpp
ParallelSolverAdapter<NuMVC> solver(...);
```

The interface of the adapter itself fulfills the LibMVC interface as well.
If not set, the adapter uses as many parallel solvers as the system provides
CPU cores.

## Compiling

All solvers are implemented as header only C++14 libraries.
Just incluce `<solver>/<solver>.hpp` in your project.
Each solver is also shipped with a standalone version. To use it call
make in the solvers folder.

For building the benchmarks, see the corresponding section below.

### Documentation

To build the documentation for all headers using doxygen, call
`make doc` in the build directory.

## How to use

The standalone version calculates the MVC of a graph given
as a text file in dimacs format.

Example using the NuMVC solver:

```bash
make
# read file frb45-21-1.mis
# stop calculation only if time is up
# use 10 seconds for calculation
./numvc frb45-21-1.mis 0 10
```

**Integration with BGL**

For integration with the Boost Graph Library (BGL) a converter is provided in
[this gist](https://gist.github.com/fmoessbauer/163b9928ae9170cfe2651173f416314b).

```cpp
// convert graph to dimacs
std::stringstream os;
write_dimacs(col_graph, os);

// initialize solver and solve
NuMVC solver(os, boost::num_vertices(col_graph)-size_is, max_sec, verbose);
solver.cover_LS();
std::vector<int> solution = std::move(solver.get_independent_set());
```

## Benchmarks

Benchmarks are provided for all solvers, using the google benchmark library.
The library is included as a git submodule, hence clone using `git clone --recursive`.
To simplify building and also as a good starting point, the `build.sh` script can be used.

To keep the repository small, only a small sample graph is included in `bench/data`.
More and larger graphs can be downloaded using `bench/fetch_graphs.sh`.

For running the benchmarks, see this example:

```bash
# Usage: LibMVC-bench <graph> <minimum-cover> <timeout (sec)> <gbench parameters>
./bench/LibMVC-bench data/frb45-21-mis/frb45-21-1.mis 900 100 --benchmark_repetitions=5
```

## Indexing of Vertices

The vertices of the DIMACS files are indexed starting at 1, however internally
the the indexing starts at 0. For convenience the methods for getting the solution
as indepentent set or cover provide a parameter to specify the indexing (default 1).

## References

For more standalone solvers see also http://lcs.ios.ac.cn/~caisw/MVC.html

