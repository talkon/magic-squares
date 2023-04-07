# magic-squares

Directory structure:
```
magic-squares
|-- research/ -> math and data analysis
|-- tests/ -> benchmarks for arrangement code
|-- data/ -> data collected from searches
|-- src/ -> 
    |-- py/ -> python code
    |-- low-level/ -> future low-level code
```

## usage

Release: `./build.sh`. Debug: `./build.sh d`.

After building, you can do `./run.sh` to run the whole pipeline.

Individual usage is in README.md of each directory.

---

Requires CMake 3.13+

Build using CMake, using `build/` as the build directory.
This generates the binary in `bin/arrangement`
```
cmake --build build
```

Build with debug info (i.e. `-g` compiler flag), useful for `perf` and the like:
```
cmake --build build --config RelWithDebInfo
```

Clean build directory:
```
cmake --build build --target clean
```
