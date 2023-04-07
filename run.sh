# usage: ./run.sh
# TODO: add real arguments...
# idk how to parse arguments in bash lmao

PRODUCT="10 4 3 2"
ARRANGEMENT_OPTIONS="--min-sum 320 --max-sum 330"
PERF=true

if [ ! -f bin/enumeration.py ]; then
  echo "run ./build.sh first!"
  exit
fi

# make enumeration file
ENUM_FILE_NAME="data/enumerations_${PRODUCT// /_}.txt"
if [ -e $ENUM_FILE_NAME ]; then
  echo "file ${ENUM_FILE_NAME} exists, not remaking it..."
else
  echo "making ${ENUM_FILE_NAME}..."
  pypy3 bin/enumeration.py $PRODUCT --file $ENUM_FILE_NAME
fi

# run arrangement
if [ "$PERF" = true ]; then
  time perf record bin/arrangement --file $ENUM_FILE_NAME $ARRANGEMENT_OPTIONS
else
  bin/arrangement --file $ENUM_FILE_NAME $ARRANGEMENT_OPTIONS
fi
