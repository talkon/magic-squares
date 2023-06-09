# usage:
#   ./build.sh [d|D] [c] [t|T]
#   release build by default
#   d: debug build (-O3)
#   D: debug build (-O0)
#   c: clean after building
#   t: test after building, no enum tests (requires enum tests to have been run once for enum files)
#   T: test after building, includes enum tests
# Note that the first run of the test cases might be slow, as enumeration files will be generated

if [[ $* == *D* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
  : ${BUILD_TYPE:=Debug}
#  cmake -D CMAKE_BUILD_TYPE= -S . -B $CMAKE_DIR
elif [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-relwithdebinfo}
  : ${BUILD_TYPE:=RelWithDebInfo}
#  cmake -D CMAKE_BUILD_TYPE=RelWithDebInfo -S . -B $CMAKE_DIR
else
  : ${CMAKE_DIR:=cmake-build-release}
  : ${BUILD_TYPE:=Release}
#  cmake -D CMAKE_BUILD_TYPE=Release -S . -B $CMAKE_DIR
fi

if [[ $* == *c* ]]; then
  rm -r $CMAKE_DIR
  rm -r bin
else
  mkdir -p $CMAKE_DIR
  cd $CMAKE_DIR
  cmake -D CMAKE_BUILD_TYPE=$BUILD_TYPE ..
  make
  cd ..
  mkdir -p bin
  ln -sf $(realpath $CMAKE_DIR)/src/c/arrangement_5 bin/arrangement_5
  ln -sf $(realpath $CMAKE_DIR)/src/c/arrangement_6 bin/arrangement_6
  ln -sf $(realpath $CMAKE_DIR)/src/c/arrangement_7 bin/arrangement_7
  ln -sf $(realpath $CMAKE_DIR)/src/c/enumeration bin/enumeration
  ln -sf $(realpath src/py/enumeration.py) bin/enumeration.py
  ln -sf $(realpath src/py/postprocess.py) bin/postprocess.py
  mkdir -p data
  mkdir -p data/5
  mkdir -p data/5/enumerations
  mkdir -p data/5/output
  mkdir -p data/5/stats
  mkdir -p data/6
  mkdir -p data/6/enumerations
  mkdir -p data/6/output
  mkdir -p data/6/stats
  mkdir -p data/7
  mkdir -p data/7/enumerations
  mkdir -p data/7/output
  mkdir -p data/7/stats
fi

if [[ $* == *t* ]]; then
  mkdir -p tests
  ctest --test-dir $CMAKE_DIR -R P
fi

if [[ $* == *T* ]]; then
  mkdir -p tests
  ctest --test-dir $CMAKE_DIR
fi
