# usage:
#   ./build.sh [d|D] [c] [t]
#   release build by default
#   d: debug build (-O3)
#   D: debug build (-O0)
#   c: clean after building
#   t: test after building
# Note that the first run of the test cases might be slow, as enumeration files will be generated

if [[ $* == *D* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
  cmake -D CMAKE_BUILD_TYPE=Debug -S . -B $CMAKE_DIR
elif [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-relwithdebinfo}
  cmake -D CMAKE_BUILD_TYPE=RelWithDebInfo -S . -B $CMAKE_DIR
else
  : ${CMAKE_DIR:=cmake-build-release}
  cmake -D CMAKE_BUILD_TYPE=Release -S . -B $CMAKE_DIR
fi

if [[ $* == *c* ]]; then
  cmake --build $CMAKE_DIR --target clean
else
  cmake --build $CMAKE_DIR
fi

mkdir -p bin
ln -sf $(realpath $CMAKE_DIR)/src/c/arrangement bin/arrangement
ln -sf $(realpath src/py/enumeration.py) bin/enumeration.py
ln -sf $(realpath src/py/postprocess.py) bin/postprocess.py

if [[ $* == *t* ]]; then
  mkdir tests
  ctest --test-dir $CMAKE_DIR
fi
