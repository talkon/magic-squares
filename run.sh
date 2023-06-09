# usage: ./run.sh 10 4 3 2 [-e[f]] [-a[f]] [-p[f]] [--perf] [--efile fname] [--afile fname] [--pfile fname] [arrangement.c options] [postprocess.py options]
# runs enumeration, arrangement, and post-processing for a given P
# options:
#   -e: run "e"numeration step
#   -a: run "a"rrangement step
#   -p: run "p"ostprocessing step
#   --glob [pattern]: run aggregated postprocessing step on all files matching specified pattern
#                     use like this: ./run.sh --glob "data/output/output_*" --pfile stats/stats_long.txt --verbose 3 -pf
#   -ef, -af, -pf: same as -e, -a, -p, but force rerunning if cached files exist
#   -c: "c"ontinue previous run
#   --perf: run perf for arrangement step
#   --vec-size: magic square dimensions, default 6
#   --output-dir [dir]: override output directory (default is data/)
#   --efile [fname]: override enumeration file name (default is data/{vec_size}/enumerations_[P].txt)
#   --afile [fname]: override arrangement file name (default is data/{vec_size}/output_[P].txt 
#                    or data/{vec_size}/output_[P]_[arrangement.c options].txt if arrangement options are used)
#   --pfile [fname]: override postprocess file name (default is data/{vec_size}/summary_[P].txt 
#                    or data/{vec_size}/summary_[P]_[arrangement.c options].txt if arrangement options are used)
# when no -e/a/p flags are used, defaults to -e -a -p
# supported arrangement.c options: --min-sum, --max-sum, --sum, --count
# supported postprocess.py options: --verbose (default 1), --cutoff

POSITIONAL_ARGS=()
ARRANGE_OPTIONS=()
POSTPROC_OPTIONS=()
MIN_SUM=0

ENUM=false
ARRANGE=false
POSTPROC=false

DEFAULT=true
PARTIAL=false

FORCE_ENUM=false
FORCE_ARRANGE=false
FORCE_POSTPROC=false

CONTINUE_ARRANGE=false
VEC_SIZE=6

GLOB=false

PERF=false
OUTPUT_DIR="data"

# allows setting environment variable for pypy3 installation
: "${PYPY3:=pypy3}"

# parse arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --min-sum)
      MIN_SUM=$2
      shift # past argument
      shift # past value
      ;;
    --max-sum|--sum|--count)
      ARRANGE_OPTIONS+=("$1" "$2")
      shift # past argument
      shift # past value
      ;;
    --expect-nsols|--cutoff)
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
    --vec-size)
      VEC_SIZE="$2"
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
    -c)
      CONTINUE_ARRANGE=true
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

if [ "$DEFAULT" = true ]; then
  ENUM=true
  ARRANGE=true
  POSTPROC=true
fi

PRODUCT="${POSITIONAL_ARGS[@]}"
ARRANGE_PATTERN="${OUTPUT_DIR}/${VEC_SIZE}/output/output_${PRODUCT// /_}.*"

if [ "$CONTINUE_ARRANGE" = true ]; then
  echo "continuing arrangement"
  # apologies for this bash one-liner
  MAX_PREV_SUM=$(cat $ARRANGE_PATTERN | grep sum | awk '{print $2}' | sort -n | tail -1)
  if [ ! "$MAX_PREV_SUM" = "" ]; then
    MIN_SUM=$(($MAX_PREV_SUM + 1))
  fi
  echo "starting at sum ${MIN_SUM}"
fi

FNAME_SUFFIX=".${MIN_SUM}"
if [ ! "$MIN_SUM" = "0" ]; then 
  ARRANGE_OPTIONS+=("--min-sum" "$MIN_SUM")
fi

ARRANGE_OPTIONS="${ARRANGE_OPTIONS[@]}"
POSTPROC_OPTIONS="${POSTPROC_OPTIONS[@]}"

if [ ! -f bin/enumeration.py ]; then
  echo "run ./build.sh first!"
  exit
fi

: ${ENUM_FILE_NAME:="${OUTPUT_DIR}/${VEC_SIZE}/enumerations/enumerations_${PRODUCT// /_}.txt"}
: ${ARRANGE_FILE_NAME:="${OUTPUT_DIR}/${VEC_SIZE}/output/output_${PRODUCT// /_}${FNAME_SUFFIX}.txt"}
: ${POSTPROC_FILE_NAME:="${OUTPUT_DIR}/${VEC_SIZE}/stats/summary_${PRODUCT// /_}${FNAME_SUFFIX}.txt"}

# make enumeration file
if [ "$ENUM" = true ]; then
  echo "running enumeration"
  if [ "$FORCE_ENUM" = false ] && [ -e $ENUM_FILE_NAME ]; then
    echo "file ${ENUM_FILE_NAME} exists, not remaking it..."
  else
    echo "making ${ENUM_FILE_NAME}..."
    $PYPY3 bin/enumeration.py $PRODUCT --file $ENUM_FILE_NAME --vec-size $VEC_SIZE
  fi
fi

# run arrangement
if [ "$ARRANGE" = true ]; then
  echo "running arrangement"
  if [ "$PERF" = true ]; then
    RUN_COMMAND="time perf record bin/arrangement_$VEC_SIZE --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
  else
    RUN_COMMAND="bin/arrangement_$VEC_SIZE --file $ENUM_FILE_NAME $ARRANGE_OPTIONS"
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
      $PYPY3 bin/postprocess.py --vec-size "$VEC_SIZE" --glob "$GLOB_STR" $POSTPROC_OPTIONS > $POSTPROC_FILE_NAME
    else
      $PYPY3 bin/postprocess.py --vec-size "$VEC_SIZE" --glob "$ARRANGE_PATTERN" --sort "P" $POSTPROC_OPTIONS > $POSTPROC_FILE_NAME && cat $POSTPROC_FILE_NAME
    fi
  fi
fi
