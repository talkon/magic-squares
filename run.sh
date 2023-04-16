# usage: ./run.sh 10 4 3 2 [-e[f]] [-a[f]] [-p[f]] [--perf] [--efile fname] [--afile fname] [--pfile fname] [arrangement.c options] [postprocess.py options]
# runs enumeration, arrangement, and post-processing for a given P
# options:
#   -e: run "e"numeration step
#   -a: run "a"rrangement step
#   -p: run "p"ostprocessing step
#   --glob [pattern]: run aggregated postprocessing step on all files matching specified pattern
#                     use like this: 
#   -ef, -af, -pf: same as -e, -a, -p, but force rerunning if cached files exist
#   --perf: run perf for arrangement step
#   --output-dir [dir]: override output directory (default is data/)
#   --efile [fname]: override enumeration file name (default is data/enumerations_[P].txt)
#   --afile [fname]: override arrangement file name (default is data/output_[P].txt or data/output_[P]_[arrangement.c options].txt if arrangement options are used)
#   --pfile [fname]: override postprocess file name (default is data/summary_[P].txt or data/summary_[P]_[arrangement.c options].txt if arrangement options are used)
# when no -e/a/p flags are used, defaults to -e -a -p
# supported arrangement.c options: --min-sum, --max-sum, --sum
# supported postprocess.py options: --verbose (default 1), --sort

POSITIONAL_ARGS=()
ARRANGE_OPTIONS=()
POSTPROC_OPTIONS=()
FNAME_SUFFIX=()

ENUM=false
ARRANGE=false
POSTPROC=false

DEFAULT=true
PARTIAL=false

FORCE_ENUM=false
FORCE_ARRANGE=false
FORCE_POSTPROC=false

GLOB=false

PERF=false
OUTPUT_DIR="data"

# allows setting environment variable for pypy3 installation
: "${PYPY3:=pypy3}"

while [[ $# -gt 0 ]]; do
  case $1 in
    --min-sum|--max-sum|--sum)
      ARRANGE_OPTIONS+=("$1" "$2")
      FNAME_SUFFIX+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --expect-nsols)
      POSTPROC_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --efile)
      ENUM_FILE_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    --afile)
      ARRANGE_FILE_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    --pfile)
      POSTPROC_FILE_NAME="$2"
      shift # past argument
      shift # past value
      ;;
    --output-dir)
      OUTPUT_DIR="$2"
      shift
      shift
      ;;
    --perf)
      PERF=true
      shift
      ;;
    --verbose)
      POSTPROC_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --sort)
      POSTPROC_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --glob)
      GLOB=true
      POSTPROC=true
      DEFAULT=false
      GLOB_STR="$2"
      shift # past argument
      shift # past value
      ;;
    -e)
      ENUM=true
      DEFAULT=false
      shift # past argument
      ;;
    -a)
      ARRANGE=true
      DEFAULT=false
      shift # past argument
      ;;
    -p)
      POSTPROC=true
      DEFAULT=false
      shift # past argument
      ;;
    -ef)
      ENUM=true
      FORCE_ENUM=true
      DEFAULT=false
      shift # past argument
      ;;
    -af)
      ARRANGE=true
      FORCE_ARRANGE=true
      DEFAULT=false
      shift # past argument
      ;;
    -pf)
      POSTPROC=true
      FORCE_POSTPROC=true
      DEFAULT=false
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      FNAME_SUFFIX+=("$1")
      shift # past argument
      ;;
  esac
done

if [ "$DEFAULT" = true ]; then
  ENUM=true
  ARRANGE=true
  POSTPROC=true
fi

PRODUCT="${POSITIONAL_ARGS[@]}"
ARRANGE_OPTIONS="${ARRANGE_OPTIONS[@]}"
POSTPROC_OPTIONS="${POSTPROC_OPTIONS[@]}"
FNAME_SUFFIX="${FNAME_SUFFIX[@]}"

if [ ! -f bin/enumeration.py ]; then
  echo "run ./build.sh first!"
  exit
fi

: ${ENUM_FILE_NAME:="${OUTPUT_DIR}/enumerations_${PRODUCT// /_}.txt"}
: ${ARRANGE_FILE_NAME:="${OUTPUT_DIR}/output_${FNAME_SUFFIX// /_}.txt"}
: ${POSTPROC_FILE_NAME:="${OUTPUT_DIR}/summary_${FNAME_SUFFIX// /_}.txt"}

# make enumeration file
if [ "$ENUM" = true ]; then
  echo "running enumeration"
  if [ "$FORCE_ENUM" = false ] && [ -e $ENUM_FILE_NAME ]; then
    echo "file ${ENUM_FILE_NAME} exists, not remaking it..."
  else
    echo "making ${ENUM_FILE_NAME}..."
    $PYPY3 bin/enumeration.py $PRODUCT --file $ENUM_FILE_NAME
  fi
fi

# run arrangement
if [ "$ARRANGE" = true ]; then
  echo "running arrangement"
  if [ "$PERF" = true ]; then
    RUN_COMMAND="time perf record bin/arrangement --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
  else
    RUN_COMMAND="bin/arrangement --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
  fi
  if [ "$FORCE_ARRANGE" = false ] && [ -e $ARRANGE_FILE_NAME ]; then
    echo "file ${ARRANGE_FILE_NAME} exists, not remaking it..."
  else
    echo "making ${ARRANGE_FILE_NAME}..."
    $RUN_COMMAND > $ARRANGE_FILE_NAME
  fi
fi

# process output
if [ "$POSTPROC" = true ]; then
  echo "running postprocessing"
  if [ "$FORCE_POSTPROC" = false ] && [ -e $POSTPROC_FILE_NAME ]; then
    echo "file ${POSTPROC_FILE_NAME} exists, not remaking it..."
    cat $POSTPROC_FILE_NAME
  else
    echo "making ${POSTPROC_FILE_NAME}..."
    if [ "$GLOB" = true ]; then
      $PYPY3 bin/postprocess.py --glob "$GLOB_STR" $POSTPROC_OPTIONS > $POSTPROC_FILE_NAME && cat $POSTPROC_FILE_NAME
    else
      $PYPY3 bin/postprocess.py --file $ARRANGE_FILE_NAME $POSTPROC_OPTIONS > $POSTPROC_FILE_NAME && cat $POSTPROC_FILE_NAME
    fi
  fi
fi
