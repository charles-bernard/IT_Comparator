#!/usr/bin/awk

BEGIN {
	FS = "\t";
	header = "Start" \
		"\tEnd" \
		"\tInitially Detect Strand" \
		"\tPrimary Sequence" \
		"\tSecondary Structure (Dot-Bracket Notation)" \
		"\tFree Energy (kcal/mol)"
	printf("%s\n", header);
}

NR > 1 {
	printf("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $6, $21, $20);
}

