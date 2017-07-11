#!/bin/bash
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
TMP_TOOL_STDERR=$(mktemp);

###########################################################
# WELCOME MESSAGE
###########################################################
cat <<WELCOME_MESSAGE
###########################################################
You are using IT Comparator!

Let us compare the genomic context of homologuous ITs 
between related bacterial species!

WELCOME_MESSAGE

###########################################################
# I. FUNCTIONS
###########################################################
###### I.1 Handle Errors ##################################
###########################################################
function error_exit {
	# exit with error message and error code sent as args
	local ERR_MSG="$1"; local ERR_CODE="$2"; local LOG="${3:-$(mktemp)}";
	printf "ERROR:\n""$ERR_MSG""\n" | tee -a "$LOG" >&2; 
	exit "${ERR_CODE:-1}";				
}

function check_tool_stderr {
	# check if a tool has printed in stderr, if yes exit.
	local TOOL_EXEC_STDERR_FILE="$1";
	local TOOL_NAME="$2"; local LOG="$3";

	if [ -s "$TOOL_EXEC_STDERR_FILE" ]; then
		local TOOL_STDERR="The tool "\""$TOOL_NAME"\"" exited with the following error:\n" 
		TOOL_STDERR="$TOOL_STDERR"`cat "$TOOL_EXEC_STDERR_FILE"`
		rm "$TOOL_EXEC_STDERR_FILE";
		error_exit "$TOOL_STDERR" 4 "$LOG";
	fi
}

###########################################################
###### I.2 Test Pipeline Parameters #######################
###########################################################
function check_param {
	local PARAM="$1"; local NAME="$2"; local OPTION="$3";
	if [ ! "$PARAM" ]; then
		ERROR="Parameter "\""$NAME"\"" is required!\n";
		ERROR="$ERROR""Use option ""$OPTION"
		error_exit "$ERROR" 1;
	fi
}

function check_indir {
	# Check conformity of Input Dir
	local DIR="$1";
	if [ ! -d "$DIR" ]; then
		error_exit \
		"The input dir \"$DIR\" does not exist!"
	fi
	if ! ls -d "$DIR"/*/ 1>$(mktemp); then
		error_exit \
		"The input dir \"$DIR\" has no sub dir.\nOne sub dir is required per species!" 3; 
	fi
}

function check_insubdir {
	# Check whether current subdir stores the expected files
	local DIR="$1";

ERR_MSG="The structure of the input sub dir \"$DIR\" does not\n"\
"satisfied the expected requirements...\n\n"\
"Any subdir must contain the three following files:\n"\
"\t* genome_sequence.(fa/fasta)\n"\
"\t* genome_annotation.(gtf/gff/gff3)\n"\
"\t* list_ITs.gff"

	if [ ! `ls "$DIR" | egrep '^genome_sequence.fa(sta)?$'` ]; then
		error_exit "$ERR_MSG";
	fi
	if [ ! `ls "$DIR" | egrep '^genome_annotation.g[tf]f(3)?$'` ]; then
		error_exit "$ERR_MSG";
	fi
	if [ ! `ls "$DIR" | egrep '^list_ITs.gff$'` ]; then
		error_exit "$ERR_MSG";
	fi
}

function check_outdir {
	# Check if output dir exists; create it if needed 
	local DIR="$1";
	if [ ! -d "$DIR" ]; then
		printf "The ouput dir \"$DIR\" does not exist.\n"; 
		printf "Do you want to create it? (y|n)\n";
		while true; do
			read -p "> " YN
			case $YN in
				[Yy]* )	mkdir -p "$DIR"; break;;
				[Nn]* ) error_exit \
					"The output dir \"$DIR\" was not valid!" 3; break;;
				* )		printf "Please answer yes or no.\n";;
			esac
		done
		printf "\n";
	fi
}

###########################################################
# II. PARAMETERS FOR THE PIPELINE
###########################################################
###### II.1 Read Options ##################################
###########################################################
TEMP=`getopt \
	-o \
		i:o:l:: \
	--long \
		input_dir:,output-dir:,log:: \
	-n 'report' -- "$@"`;
eval set -- "$TEMP";

###########################################################
###### II.2 Extract Options ###############################
###########################################################
while true ; do
	case "$1" in
		-i | --input-dir )
			case "$2" in
				"") error_exit "No input dir has not been provided!" 1; shift 2 ;;
				*) INPUT_DIR="$2"; shift 2 ;;
			esac ;;
		-o | --output-dir )
			case "$2" in
				"") error_exit "No output dir has not been provided!" 1; shift 2 ;;
				*) OUTPUT_DIR="$2"; shift 2 ;;
			esac ;;
		-l | --log )
			case "$2" in
				*) LOG="$2"; shift 2 ;;
			esac ;;
		-- ) shift; break ;;
		*) error_exit "Internal error!" 1; 
	esac
done

###########################################################
###### II.3 Check pipeline parameters #####################
###########################################################
printf "###########################################################\n"
printf "PRE-PROCESSING)\n";
printf " * Checking the validity of the pipeline parameters\n\n";

RECQ_PARAMS=("$INPUT_DIR" "$OUTPUT_DIR");
PARAM_NAMES=("Input Directory" "Output Directory");
OPTIONS=("-i/--input-dir" "-o/--output-dir");
for (( i=0; i<${#RECQ_PARAMS[@]}; i++ )); do 
	check_param "${RECQ_PARAMS[$i]}" "${PARAM_NAMES[$i]}" "${OPTIONS[$i]}";
done

check_indir "$INPUT_DIR";
s=0;
for i in $(ls -d "$INPUT_DIR"/*/); do
	INSUB_DIR[$s]="$(basename $i)"; 
	check_insubdir "$INPUT_DIR"/"${INSUB_DIR[$s]}"
	s=$(($s+1));
done


check_outdir "$OUTPUT_DIR"; 

printf "Parameters Check was successful!\n\n";

COMMAND="$OUTPUT_DIR"/"commands.log"; > "$COMMAND";