#!/usr/bin/awk

BEGIN {
	FS = "\t";
	file_idx = 0;
	# aligner is sent to this script as arg
}

FNR == 1 {
	file_idx++;
}

file_idx == 1 { 
	ref_seq[$1] = $2;
}

file_idx == 2 {
	vs_seq[$1] = $2;
}

file_idx == 3 {
	aligned = $19;

	if(aligned == "") {
		alignment = "";
		if(vs_seq[$3]) {
			cmd = aligner " " ref_seq[$1] " " vs_seq[$3];
			cmd |& getline alignment
		}

		# print updated line
		for(i = 1; i < 19; i++) {
			printf("%s\t", $i);
		}
		printf("%s\t%s\t%s\n", alignment, $20, $21);

	} else {
		printf("%s\n", $0);
	}
}