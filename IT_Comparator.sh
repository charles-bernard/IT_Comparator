#!/bin/bash
SCRIPT_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )";
TMP_TOOL_STDERR=$(mktemp);

###########################################################
# WELCOME MESSAGE
###########################################################
cat <<WELCOME_MESSAGE
###########################################################
You are using IT Comparator!

Let us compare the genomic landscape/context of 
homologuous ITs among related bacterial species!

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
	local DIR="$1"; local ISREF="${2:-0}"
	if [ ! -d "$DIR" ]; then
		error_exit \
		"The input dir \"$DIR\" does not exist!"
	fi
	if [ "$ISREF" == 0 ]; then
		if ! ls -d "$DIR"/*/ 1>$(mktemp); then
			error_exit \
			"The input dir \"$DIR\" has no sub dir.\nOne sub dir is required per species!" 3; 
		fi
	fi
}

function check_fasta {
	local GENOME="$1"; local DIR="$2";
	read -r -n 1 FIRST_CHAR < "$DIR"/"$GENOME";
	if [[ "$FIRST_CHAR" != ">" ]]; then
		error_exit \
		"The fasta file in \"$DIR\" does\nnot begin with the required \">\" character!" 2;
	fi
}

function check_insubdir {
	# Check whether current subdir stores the expected files
	local DIR="$1";

ERR_MSG="The structure of the input sub dir \"$DIR\"\n"\
"does not satisfy the expected requirements...\n\n"\
"Any subdir must contain the three following files:\n"\
"\t* genome_sequence.(fa/fasta)\n"\
"\t* genome_annotation.(gtf/gff/gff3)\n"\
"\t* list_ITs.gff"

	FASTA=`ls "$DIR" | egrep '^genome_sequence.fa(sta)?$'`
	if [ ! "$FASTA" ]; then
		error_exit "$ERR_MSG";
	fi
	check_fasta "$FASTA" "$DIR";
	if [ ! `ls "$DIR" | egrep '^genome_annotation.g[tf]f(3)?$'` ]; then
		error_exit "$ERR_MSG";
	fi
	if [ ! `ls "$DIR" | egrep '^list_ITs.gff$'` ]; then
		error_exit "$ERR_MSG";
	fi
}

function confirm_species_name {
	CUR_SPECIES_NAME="$1";

	printf "By default, the name of the species is given\n";
	printf "by the name of their directory\n";
	while true; do
		printf "Do you confirm \"$CUR_SPECIES_NAME\" as species name? (y|n)\n"
		read -p "> " YN
			case $YN in
				[Yy]* )	break;;
				[Nn]* ) printf "Please, type the name of the species\n"; \
						read -p "> " CUR_SPECIES_NAME; break;;
				* )		printf "Please answer yes or no.\n";;
			esac
	done
	printf "\"$CUR_SPECIES_NAME\" saved successfully as species name!\n\n";
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

function check_log {
	# Create log file if necessary
	local LOCAL_LOG="$1";
	if [ ! $LOCAL_LOG ]; then 
		LOG="$OUTPUT_DIR"/IT_Comparator.log;
		printf "You did not mention a path for the log file.\n";
		printf "The log file will be automatically generated at:\n";
		printf "\t\"$LOG\"\n\n";
	fi
}

###########################################################
###### I.3 GFF to tab #####################################
###########################################################
function gff_to_tab {
	# convert output of IT_Miner to csv

	# Args
	local SCRIPT="$1"; local IN_FILE="$2";
	local OUT_FILE="$3"; local COMMAND="$4";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "CONVERT IT_Miner GFF OUTPUT TO TABULAR FILE:\n" >> "$COMMAND";

	(set -x; 
		awk -f "$SCRIPT" "$IN_FILE" \
		> "$OUT_FILE" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.4 translate ref genes ############################
###########################################################
function translate_ref_genes {
	# get protein seq of ref genes

	# Args
	local SCRIPT="$1"; local GENOME="$2";
	local ANNOTATION="$3"; OUT_FILE="$4"; local COMMAND="$5";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "TRANSLATE REF GENES:\n" >> "$COMMAND";

	(set -x; 
		awk -f "$SCRIPT" "$GENOME" "$ANNOTATION" \
		> "$OUT_FILE" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
# II. PARAMETERS FOR THE PIPELINE
###########################################################
###### II.1 Read Options ##################################
###########################################################
TEMP=`getopt \
	-o \
		i:r:o:l:: \
	--long \
		input_dir:,reference-dir:,output-dir:,log:: \
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
		-r | --reference-dir )
			case "$2" in
				"") REF_DIR="$SCRIPT_PATH"/"Default_reference"/"Escherichia_coli_K12_U00096.2"; shift 2 ;;
				*) REF_DIR="$2"; shift 2 ;;
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

if [ ! "$REF_DIR" ]; then 
	REF_DIR="$SCRIPT_PATH"/"Default_reference"/"Escherichia_coli_K12_U00096.2"; 
fi
check_indir "$REF_DIR" 1; check_insubdir "$REF_DIR";
confirm_species_name "$(basename "$REF_DIR")"; REF_SPECIES=$CUR_SPECIES_NAME;

check_indir "$INPUT_DIR";
s=0;
for i in $(ls -d "$INPUT_DIR"/*/); do
	INSUB_DIR[$s]="$(basename $i)"; 
	check_insubdir "$INPUT_DIR"/"${INSUB_DIR[$s]}";
	confirm_species_name "${INSUB_DIR[$s]}";
	SPECIES_NAME[$s]=$CUR_SPECIES_NAME;
	if [ "${INSUB_DIR[$s]}" != "$REF_DIR" ]; then
		s=$(($s+1));
	fi
done
NB_SPECIES=$s;

check_outdir "$OUTPUT_DIR"; 

check_log "$LOG"; > "$LOG";
COMMAND="$OUTPUT_DIR"/"commands.log"; > "$COMMAND";

printf "Parameters Check was successful!\n\n";

###########################################################
###### II.4 Print pipeline parameters #####################
###########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "PARAMETERS LIST)\n" | tee -a "$LOG";
printf " * OUTPUT DIRECTORY:"\ \""$OUTPUT_DIR"\""\n" | tee -a "$LOG";
printf " * REF DIRECTORY:   "\ \""$REF_DIR"\""\n" | tee -a "$LOG";
printf " * INPUT DIRECTORY: "\ \""$INPUT_DIR"\""\n" | tee -a "$LOG";
printf " * LOG FILE:        "\ \""$LOG"\""\n" | tee -a "$LOG";
printf " * REF SPECIES:     "\ \""$REF_SPECIES"\""\n" | tee -a "$LOG";
printf " * NB_SPECIES:      "" $NB_SPECIES""\n" | tee -a "$LOG";
for (( i=0; i<$NB_SPECIES; i++ )); do
	printf "   * SPECIES $i:     "\
\ \""${SPECIES_NAME[$i]}"\""\n" | tee -a "$LOG";
done

##########################################################
##### III. Prepare Input #################################
##########################################################
###### III.1 Create outdirs ##############################
##########################################################
REF_OUTDIR="$OUTPUT_DIR"/"Reference";
mkdir -p "$REF_OUTDIR";
for (( i=0; i<${#SPECIES_NAME[@]}; i++ )) do
	SPECIES_OUTDIR[$i]="$OUTPUT_DIR"/"${SPECIES_NAME[$i]}";
	mkdir -p "${SPECIES_OUTDIR[$i]}";
done
##########################################################
###### GFF To Tab ########################################
##########################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "PREPARING INPUT DATA\n" | tee -a "$LOG";
printf "STEP 00) list_ITs.gff to tabular file\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"00-gff_to_tab.awk";

REF_IT_GFF="$REF_DIR"/"list_ITs.gff";
REF_TAB_FILE="$REF_OUTDIR"/"00-list_ITs.csv";
gff_to_tab "$SCRIPT" "$REF_IT_GFF" "$REF_TAB_FILE" "$COMMAND";

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	IT_GFF="$INPUT_DIR"/"${INSUB_DIR[$i]}"/"list_ITs.gff";
	TAB_FILE[$i]="${SPECIES_OUTDIR[$i]}"/"00-list_ITs.csv";
	gff_to_tab "$SCRIPT" "$IT_GFF" "${TAB_FILE[$i]}" "$COMMAND";
done
##########################################################
###### Translate Ref genes ###############################
##########################################################
printf "STEP 01) Translate Reference Genes\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"01-translate_ref_genes.awk";

REF_GENOME_SEQ="$REF_DIR"/"`ls "$REF_DIR" | egrep '^genome_sequence.fa(sta)?$'`";
REF_GENOME_ANNOTATIONS="$REF_DIR"/"`ls "$REF_DIR" | egrep '^genome_annotation.g[tf]f(3)?$'`";
REF_GENES="$REF_OUTDIR"/"01-ref_genes.csv";

translate_ref_genes "$SCRIPT" "$REF_GENOME_SEQ" \
	"$REF_GENOME_ANNOTATIONS" "$REF_GENES" "$COMMAND"; 