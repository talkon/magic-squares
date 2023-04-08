# usage:
#   ./build.sh: build release version
#   ./build.sh d: build debug version
#   ./build.sh c: clean release build
#   ./build.sh cd: clean debug build

if [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
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
ln -sf $(realpath src/py/verifier.py) bin/verifier.py
