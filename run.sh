# usage: ./run.sh 10 4 3 2 [-p] [-e | -a | -l] [--output] [arrangement.c options] [postprocess.py options]
# runs enumeration, arrangement, and post-processing for a given P
# options:
#   -e: only run "e"numeration step
#   -a: run up to "a"rrangement step (i.e. no post-processing)
#   -l: "l"azy, only do post-processing step if arrangement output already exists
#   -p: run "p"erf for arrangement step
#   --output: override output file name
# supported arrangement.c options: --min-sum, --max-sum, --sum
# supported postprocess.py options: --expect-nsols

POSITIONAL_ARGS=()
ARRANGE_OPTIONS=()
POSTPROC_OPTIONS=()

ENUM=true
FORCE_ENUM=false
ARRANGE=true
FORCE_ARRANGE=true
POSTPROC=true

PERF=false

# allows setting environment variable for pypy3 installation
: "${PYPY3:=pypy3}"

while [[ $# -gt 0 ]]; do
  case $1 in
    --min-sum|--max-sum|--sum)
      ARRANGE_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --expect-nsols)
      POSTPROC_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --output)
      OUTPUT_FILE_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    -p)
      PERF=true
      shift # past argument
      ;;
    -e)
      ARRANGE=false
      POSTPROC=false
      shift # past argument
      ;;
    -a)
      POSTPROC=false
      shift # past argument
      ;;
    -l)
      FORCE_ARRANGE=false
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
ARRANGE_OPTIONS="${ARRANGE_OPTIONS[@]}"
POSTPROC_OPTIONS="${POSTPROC_OPTIONS[@]}"

if [ ! -f bin/enumeration.py ]; then
  echo "run ./build.sh first!"
  exit
fi

# make enumeration file
if [ "$ENUM" = true ]; then
  ENUM_FILE_NAME="data/enumerations_${PRODUCT// /_}.txt"
  if [ "$FORCE_ENUM" = false ] && [ -e $ENUM_FILE_NAME ]; then
    echo "file ${ENUM_FILE_NAME} exists, not remaking it..."
  else
    echo "making ${ENUM_FILE_NAME}..."
    $PYPY3 bin/enumeration.py $PRODUCT --file $ENUM_FILE_NAME
  fi
fi

# run arrangement
if [ "$ARRANGE" = true ]; then
  if [ "$PERF" = true ]; then
    RUN_COMMAND="time perf record bin/arrangement --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
  else
    RUN_COMMAND="bin/arrangement --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
  fi
  : ${OUTPUT_FILE_NAME:="data/output_${PRODUCT// /_}.txt"}
  if [ "$FORCE_ARRANGE" = false ] && [ -e $OUTPUT_FILE_NAME ]; then
    echo "file ${OUTPUT_FILE_NAME} exists, not remaking it..."
  else
    echo "making ${OUTPUT_FILE_NAME}..."
    $RUN_COMMAND > $OUTPUT_FILE_NAME
  fi
fi

# process output
if [ "$POSTPROC" = true ]; then
  $PYPY3 bin/postprocess.py $OUTPUT_FILE_NAME -v $POSTPROC_OPTIONS
fi
