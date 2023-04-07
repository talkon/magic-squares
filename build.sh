# usage: ./build.sh or ./build.sh d

if [[ $* == *d* ]]; then
  : ${CMAKE_DIR:=cmake-build-debug}
else
  : ${CMAKE_DIR:=cmake-build-release}
fi

cmake -D CMAKE_BUILD_TYPE=Release -S . -B $CMAKE_DIR
cmake --build $CMAKE_DIR

ln -sf $(realpath $CMAKE_DIR)/src/c/arrangement bin/arrangement
ln -sf $(realpath src/py/enumeration.py) bin/enumeration.py
