#!/usr/bin/awk

BEGIN {
	FS = "\t";
	file_idx = 0;
	# var 'list_species' is sent to this script
	# var 'it_file' and 'code_file' also

	printf("Ref IT\t") > code_file;
	split(list_species, species, ";");
	k=1;
	while(species[k+1]) {
		printf("%s\t", species[k]) > it_file;
		printf("%s\t", species[k]) > code_file;
		k++;
	}
	printf("%s\n", species[k]) > it_file;
	printf("%s\n", species[k]) > code_file;
}

FNR == 1 {
	file_idx++;
	n = 1;
}

file_idx == 1 {
	ref_it[n] = $1;
	ref_code[n] = $2;
}

FNR > 1 {
	it[file_idx][n] = $3;
	if($5 && $6) {
		code[file_idx][n] = $5 $6
	} else if($5) {
		code[file_idx][n] = $5 "!"
	} else if($6) {
		code[file_idx][n] = "!" $6
	} 
	n++;
}

END  {
	for(i = 1; i <= n; i++) {
		printf("%s", ref_it[i]) > it_file;
		printf("%s\t%s", ref_it[i], ref_code[i]) > code_file ;
		for(j = 1; j <= file_idx; j++) {
			printf("\t%s", it[j][i]) > it_file;
			printf("\t%s", code[j][i]) > code_file;
		}
		printf("\n") > it_file;
		printf("\n") > code_file;
	}
}