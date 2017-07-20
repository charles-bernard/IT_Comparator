#!/usr/bin/awk

BEGIN {
	FS = "\t";

	file_idx = 0;

	# The arg 'column' is sent to this script
}

FNR == 1 {
	file_idx++;

}

file_idx == 1 && FNR > 1 && $column != "." {
	ref_name[$column] = $1;
}

file_idx == 2 && FNR == 1 {
	header = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5;
	printf("%s\n", header);
}

file_idx == 2 && FNR > 1 {
	gene = $4;
	if(ref_name[gene]) {
		printf("%s\t%s\t%s\tH:%s\t%s\n", $1, $2, $3, ref_name[gene], $5); 
	} else {
		printf("%s\t%s\t%s\t?:%s\t%s\n", $1, $2, $3, gene, $5); 
	}
}