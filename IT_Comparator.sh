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
homologous ITs among related bacterial species!

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
function translate_genes {
	# get protein seq of ref genes

	# Args
	local SCRIPT="$1"; local GENOME="$2";
	local ANNOTATION="$3"; OUT_PREFIX="$4"; local COMMAND="$5";

	local ALL_GENES="$OUT_PREFIX""all_genes.csv";
	local CDS_FILE="$OUT_PREFIX""genes-cds.csv";
	local ABN_CDS_FILE="$OUT_PREFIX""genes-abnormal_cds.csv";
	local NC_FILE="$OUT_PREFIX""genes-non_coding.csv";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "TRANSLATE REF GENES:\n" >> "$COMMAND";

	(set -x; 
		awk -v cds_file="$CDS_FILE" \
		-v abn_cds_file="$ABN_CDS_FILE" \
		-v nc_file="$NC_FILE" \
		-f "$SCRIPT" "$GENOME" "$ANNOTATION" \
		> "$ALL_GENES" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.5 Build Homology Tables ##########################
###########################################################

function homolo_build {
	# build homolo table between ref species and target species

	# Args
	local SCRIPT="$1"; local REF_FILE="$2";
	local TITLE="$3"; local VS_FILE="$4";
	local OUT_FILE="$5"; local COMMAND="$6";
	local TMP_TOOL_STDERR="$7"; local LOG="$8";
	local MODE="$9"; 
	local OPTA="${10:-""}"; local OPTB="${11:-""}"; 

	printf "\t* ""$TITLE""\n" | tee -a "$LOG";

	(set -x; 
		"$SCRIPT" -r "$REF_FILE" \
		-c "$VS_FILE" -o "$OUT_FILE" \
		-m "$MODE" $OPTA $OPTB \
		2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.6 Correct Tables #################################
###########################################################
function correct_table {
	# Args
	local SCRIPT="$1"; local INFILE="$2";
	local COMMAND="$3";
	local OUTFILE="$INFILE""_corrected";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "CORRECT TABLE:\n" >> "$COMMAND";

	(set -x; 
		awk -f "$SCRIPT" "$INFILE" \
		> "$OUTFILE" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";
	mv "$OUTFILE" "$INFILE";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.7 Modify Annotation ##############################
###########################################################
function modify_annotation {
	# Args
	local SCRIPT="$1"; local COL="$2";
	local HOMOLO_TAB="$3"; local OLD_ANNOTATION="$4"; 
	local OUTFILE="$5"; local COMMAND="$6";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "MODIFY ANNOTATION:\n" >> "$COMMAND";

	(set -x;
		awk -v column="$COL" \
		-f "$SCRIPT" \
		"$HOMOLO_TAB" "$OLD_ANNOTATION" \
		> "$OUTFILE" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.8 Get Genomic Landscape ##########################
###########################################################

function get_landscape {
	# ARGS
	local SCRIPT="$1"; BREAK_IT=$2;
	local ANNOTATION="$3"; local IT_FILE="$4"; 
	local OUTFILE="$5"; local COMMAND="$6";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "GET GENOMIC LANDSCAPE:\n" >> "$COMMAND";

	(set -x;
		awk -v breakIT=$BREAK_IT \
		-f "$SCRIPT" \
		"$ANNOTATION" "$IT_FILE" \
		> "$OUTFILE" 2>"$TMP_TOOL_STDERR";
	) 2>> "$COMMAND";
	printf "\n" >> "$COMMAND";

	# Exit if stderr_file not empty
	check_tool_stderr "$TMP_TOOL_STDERR" "$(basename "$SCRIPT")" "$LOG";
}

###########################################################
###### I.9 Get Homology of Genomic Landscape ##############
###########################################################

function get_homolo_landscape {
	# Args
	local SCRIPT="$1"; local LANDSCAPE="$2";
	local REF_LANDSCAPE="$3"; local OUTFILE="$4";
	local COMMAND="$5";

	printf "___________________________________________________________\n" >> "$COMMAND"
	printf "GET GENOMIC LANDSCAPE HOMOLOGY:\n" >> "$COMMAND";

	(set -x;
		awk -f "$SCRIPT" \
		"$LANDSCAPE" "$REF_LANDSCAPE" \
		> "$OUTFILE" 2>"$TMP_TOOL_STDERR";
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
	SPECIES_DIR[$s]="$INPUT_DIR"/"${INSUB_DIR[$s]}";
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
REF_IT_TAB_FILE="$REF_OUTDIR"/"00-list_ITs.csv";
gff_to_tab "$SCRIPT" "$REF_IT_GFF" "$REF_IT_TAB_FILE" "$COMMAND";

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	IT_GFF="${SPECIES_DIR[$i]}"/"list_ITs.gff";
	IT_TAB_FILE[$i]="${SPECIES_OUTDIR[$i]}"/"00-list_ITs.csv";
	gff_to_tab "$SCRIPT" "$IT_GFF" "${IT_TAB_FILE[$i]}" "$COMMAND";
done

##########################################################
###### Translate Genes ###################################
##########################################################
printf "STEP 01) Translate Coding Genes\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"01-translate_genes.awk";

REF_GENOME_SEQ="$REF_DIR"/"`ls "$REF_DIR" | egrep '^genome_sequence.fa(sta)?$'`";
REF_GFF_ANNOTATION="$REF_DIR"/"`ls "$REF_DIR" | egrep '^genome_annotation.g[tf]f(3)?$'`";
REF_GENES_OUTDIR="$REF_OUTDIR"/"01-Genes";
mkdir -p "$REF_GENES_OUTDIR";
REF_PREFIX_GENES="$REF_GENES_OUTDIR"/;
printf "\t* ""$REF_SPECIES""\n" | tee -a "$LOG";
translate_genes "$SCRIPT" "$REF_GENOME_SEQ" \
	"$REF_GFF_ANNOTATION" "$REF_PREFIX_GENES" "$COMMAND"; 

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	GENOME_SEQ[$i]="${SPECIES_DIR[$i]}"/"`ls "$REF_DIR" | egrep '^genome_sequence.fa(sta)?$'`";
	GFF_ANNOTATION[$i]="${SPECIES_DIR[$i]}"/"`ls "${SPECIES_DIR[$i]}" | egrep '^genome_annotation.g[tf]f(3)?$'`";
	GENES_OUTDIR[$i]="${SPECIES_OUTDIR[$i]}"/"01-Genes";
	mkdir -p "${GENES_OUTDIR[$i]}";
	PREFIX_GENES[$i]="${GENES_OUTDIR[$i]}"/;
	printf "\t* ""${SPECIES_NAME[$i]}""\n" | tee -a "$LOG";
	translate_genes "$SCRIPT" "${GENOME_SEQ[$i]}" \
		"${GFF_ANNOTATION[$i]}" "${PREFIX_GENES[$i]}" "$COMMAND"; 
done

SUFFIX_GENES[0]="genes-cds.csv";
SUFFIX_GENES[1]="genes-abnormal_cds.csv";
SUFFIX_GENES[2]="genes-non_coding.csv";
SUFFIX_ALL_GENES="all_genes.csv"

##########################################################
###### Build Homology Tables #############################
##########################################################
printf "STEP 02) Build Homology Tables of Genes\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"HomoloBuild"/"homoloBuild";
K=0; N=0;

printf "   * $REF_SPECIES _VS_ OTHER SPECIES\n";
printf "     Parallelizing the process...\n";
printf "      (Expect ~10min per Species)\n"
printf "            (What about a tea?)\n";

for (( i=0; i<${#GENES_OUTDIR[@]}; i++ )); do
	HOMOLO_OUTDIR[$i]="${SPECIES_OUTDIR[$i]}"/"02-Homology_tables";
	mkdir -p "${HOMOLO_OUTDIR[$i]}";
	for (( j=0; j<${#SUFFIX_GENES[@]}; j++ )); do
		if [ -f "${GENES_OUTDIR[$i]}"/"${SUFFIX_GENES[$j]}" ] \
		&& [ -f "$REF_GENES_OUTDIR"/"${SUFFIX_GENES[$j]}" ]; then
			HOMOLO[$N]="${HOMOLO_OUTDIR[$i]}"/"Homology-""${SUFFIX_GENES[$j]}";
			OUT_FILE[$K]="${HOMOLO[$N]}";
			REF_FILE[$K]="$REF_GENES_OUTDIR"/"${SUFFIX_GENES[$j]}";
			VS_FILE[$K]="${GENES_OUTDIR[$i]}"/"${SUFFIX_GENES[$j]}";
			TITLE[$K]="${SPECIES_NAME[$i]}"": (""${SUFFIX_GENES[$j]%.*}"")";
			K=$(($K+1));
		else
			HOMOLO[$N]="";
		fi
		N=$(($N+1));
	done
done

export -f homolo_build; export -f check_tool_stderr;
printf "___________________________________________________________\n" >> "$COMMAND"
printf "BUILD GENES HOMOLOGY:\n" >> "$COMMAND";
# parallel -k homolo_build {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} \
# 	::: "$SCRIPT" \
# 	::: "${REF_FILE[@]}" :::+ "${TITLE[@]}" :::+ "${VS_FILE[@]}" :::+ "${OUT_FILE[@]}" \
# 	::: "$COMMAND" :::+ "$TMP_TOOL_STDERR" :::+ "$LOG" :::+ "NW" :::+ "-a";

##########################################################
##### Check subtables and build reference homology tab ###
##########################################################
printf "STEP 03) Correct Tables and build Reference Table\n" | tee -a "$LOG";
SCRIPT="$SCRIPT_PATH"/"Subscripts"/"02-filter_homolotab.awk";

for(( i=0; i < K; i++ )); do
	correct_table "$SCRIPT" "${OUT_FILE[$i]}" "$COMMAND";
done

CONC_FILE[0]="$REF_PREFIX_GENES""conc_tables.csv";
printf "$REF_SPECIES""\n" > "${CONC_FILE[0]}";
for (( j=0; j<${#SUFFIX_GENES[@]}; j++ )); do
	if [ -f "$REF_PREFIX_GENES""${SUFFIX_GENES[$j]}" ]; then
		awk '{ print $1; }' \
			"$REF_PREFIX_GENES""${SUFFIX_GENES[$j]}" >> "${CONC_FILE[0]}";
	fi
done

N=0;
for (( i=0; i<${#HOMOLO_OUTDIR[@]}; i++ )); do
	CONC_FILE[$((i+1))]="${HOMOLO_OUTDIR[$i]}"/"conc_tables.csv";
	printf "${SPECIES_NAME[$i]}""\n" > "${CONC_FILE[$((i+1))]}";
	for (( j=0; j<${#SUFFIX_GENES[@]}; j++ )); do 
		if [ -f "${HOMOLO[$N]}" ]; then
			awk 'NR > 1 { print $2; }' "${HOMOLO[$N]}" >> "${CONC_FILE[$((i+1))]}";
		elif [ -f "$REF_PREFIX_GENES""${SUFFIX_GENES[$j]}" ]; then
			awk 'NR > 1 { print "."; }' \
			"$REF_PREFIX_GENES""${SUFFIX_GENES[$j]}" >> "${CONC_FILE[$((i+1))]}";
		fi
		N=$((N+1));
	done
done

mkdir -p "$OUTPUT_DIR"/"Gene_homology_table"
HOMOLO_TAB="$OUTPUT_DIR"/"Gene_homology_table"/"Homology_table.csv";
TMP=$(mktemp);
# Concatenate columns
pr -mts "${CONC_FILE[@]}" > "$TMP";
# And sort resulting table
head -1 "$TMP" > "$HOMOLO_TAB"; 
tail -n +2 "$TMP" | sort -k1 -n >> "$HOMOLO_TAB";
rm "$TMP";

##########################################################
# Modify Species Annotation according to Reference Names #
##########################################################
printf "STEP 04) Modify Annotations with respect to Reference Gene Names\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"03-modify_annotation.awk";

for ((i = 0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	OLD_TAB_ANNOTATION="${PREFIX_GENES[$i]}""$SUFFIX_ALL_GENES";
	NEW_TAB_ANNOTATION[$i]="${SPECIES_OUTDIR[$i]}"/"03-Modified_genome_annotation.csv";
	modify_annotation "$SCRIPT" "$((i+2))" "$HOMOLO_TAB" \
		"$OLD_TAB_ANNOTATION" "${NEW_TAB_ANNOTATION[$i]}" "$COMMAND";
done

##########################################################
###### Get Genomic Landscape of Genes ####################
##########################################################
printf "STEP 05) Get genomic landscape of each IT\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"04-get_genomic_landscape.awk";

REF_TAB_ANNOTATION="$REF_PREFIX_GENES""$SUFFIX_ALL_GENES";
REF_LANDSCAPE="$REF_OUTDIR"/"02-ITs_genomic_landscape.csv";
get_landscape "$SCRIPT" 1 "$REF_TAB_ANNOTATION" \
	"$REF_IT_TAB_FILE" "$REF_LANDSCAPE" "$COMMAND";

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	LANDSCAPE[$i]="${SPECIES_OUTDIR[$i]}"/"04-ITs_genomic_landscape.csv";
	get_landscape "$SCRIPT" 0 "${NEW_TAB_ANNOTATION[$i]}" \
		"${IT_TAB_FILE[$i]}" "${LANDSCAPE[$i]}" "$COMMAND";
done

##########################################################
###### Prepare ITs Seq and Struct for HomoloBuild ########
##########################################################
# TO DO: I was just lazy to integrate this part in step 00 but
# this would prevent from reading several time the same info
printf "STEP 06) Prepare list of ITs for HomologSearch\n";

REF_IT_SEQ_DIR="$REF_OUTDIR"/"03-ITs_Seq_and_Struct";
mkdir -p "$REF_IT_SEQ_DIR";
REF_IT_SEQ="$REF_IT_SEQ_DIR"/"ITs_seq.csv";
awk 'BEGIN{ FS="\t"; } NR>1{ print $1 "\t" $5; }' \
	"$REF_IT_TAB_FILE" > "$REF_IT_SEQ"; 
REF_IT_STRUCT="$REF_IT_SEQ_DIR"/"ITs_struct.csv";
awk 'BEGIN{ FS="\t"; } NR>1{ print $1 "\t" $6; }' \
	"$REF_IT_TAB_FILE" > "$REF_IT_STRUCT";

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	IT_SEQ_DIR[$i]="${SPECIES_OUTDIR[$i]}"/"05-ITs_Seq_and_Struct";
	mkdir -p "${IT_SEQ_DIR[$i]}";
	IT_SEQ[$i]="${IT_SEQ_DIR[$i]}"/"ITs_seq.csv";
	awk 'BEGIN{ FS="\t"; } NR>1{ print $1 "\t" $5; }' \
		"${IT_TAB_FILE[$i]}" > "${IT_SEQ[$i]}"; 
	IT_STRUCT[$i]="${IT_SEQ_DIR[$i]}"/"ITs_struct.csv";
	awk 'BEGIN{ FS="\t"; } NR>1{ print $1 "\t" $6; }' \
		"${IT_TAB_FILE[$i]}" > "${IT_STRUCT[$i]}"; 
done

printf "\nPreparation of Input Data was successful!\n" | tee -a "$LOG";

##############################################################
###### Get Homology of terminators by sequence & structure ###
##############################################################
printf "###########################################################\n" | tee -a "$LOG"
printf "FIND HOMOLOGOUS ITs\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"HomoloBuild"/"homoloBuild";

for (( i=0; i<${#SPECIES_OUTDIR[@]}; i++ )); do
	IT_HOMOLO_DIR[$i]="${SPECIES_OUTDIR[$i]}"/"06-ITs_Homology";
	mkdir -p "${IT_HOMOLO_DIR[$i]}";
	IT_HOMOLO_SEQ[$i]="${IT_HOMOLO_DIR[$i]}"/"Sequence_homology.csv";
	IT_HOMOLO_STRUCT[$i]="${IT_HOMOLO_DIR[$i]}"/"Structure_homology.csv";
done

printf "STEP 00) Get Sequence Homology\n" | tee -a "$LOG";
printf "___________________________________________________________\n" >> "$COMMAND"
printf "BUILD IT SEQUENCE HOMOLOGY:\n" >> "$COMMAND";
parallel -k homolo_build {1} {2} {3} {4} {5} {6} {7} {8} {9} \
	::: "$SCRIPT" :::+ "$REF_IT_SEQ" \
	::: "${SPECIES_NAME[@]}" :::+ "${IT_SEQ[@]}" :::+ "${IT_HOMOLO_SEQ[@]}" \
	::: "$COMMAND" :::+ "$TMP_TOOL_STDERR" :::+ "$LOG" :::+ "HW";

printf "STEP 01) Get Secondary Structure Homology\n" | tee -a "$LOG";
printf "___________________________________________________________\n" >> "$COMMAND"
printf "BUILD IT STRUCTURE HOMOLOGY:\n" >> "$COMMAND";
parallel -k homolo_build {1} {2} {3} {4} {5} {6} {7} {8} {9} \
	::: "$SCRIPT" :::+ "$REF_IT_STRUCT" \
	::: "${SPECIES_NAME[@]}" :::+ "${IT_STRUCT[@]}" :::+ "${IT_HOMOLO_STRUCT[@]}" \
	::: "$COMMAND" :::+ "$TMP_TOOL_STDERR" :::+ "$LOG" :::+ "HW";

##############################################################
###### Get Homology of terminators by genomic landscape ######
##############################################################
printf "STEP 02) Get Genomic Landscape Homology\n" | tee -a "$LOG";

SCRIPT="$SCRIPT_PATH"/"Subscripts"/"05-get_homologuous_landscape.awk";

for (( i=0; i<${#IT_HOMOLO_DIR[@]}; i++ )); do
	IT_HOMOLO_LANDSCAPE[$i]="${IT_HOMOLO_DIR[$i]}"/"Genomic_landscape_homology"; 
	printf "\t* ""${SPECIES_NAME[$i]}""\n" | tee -a "$LOG";
	get_homolo_landscape "$SCRIPT" "${LANDSCAPE[$i]}" \
		"$REF_LANDSCAPE" "${IT_HOMOLO_LANDSCAPE[$i]}" "$COMMAND";
done