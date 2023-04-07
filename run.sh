# usage: ./run.sh 10 4 3 2 --sum 327
# supported arrangement options: --min-sum, --max-sum, --sum

POSITIONAL_ARGS=()
ARRANGEMENT_OPTIONS=()
PERF=false

# allows setting environment variable for pypy3 installation
: "${PYPY3:=pypy3}"

while [[ $# -gt 0 ]]; do
  case $1 in
    --min-sum)
      ARRANGEMENT_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --max-sum)
      ARRANGEMENT_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --sum)
      ARRANGEMENT_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    -p)
      PERF=true
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

PRODUCT="${POSITIONAL_ARGS[@]}"
ARRANGEMENT_OPTIONS="${ARRANGEMENT_OPTIONS[@]}"

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
  $PYPY3 bin/enumeration.py $PRODUCT --file $ENUM_FILE_NAME
fi

# run arrangement
if [ "$PERF" = true ]; then
  time perf record bin/arrangement --file $ENUM_FILE_NAME $ARRANGEMENT_OPTIONS
else
  bin/arrangement --file $ENUM_FILE_NAME $ARRANGEMENT_OPTIONS
fi
