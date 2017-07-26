#!/usr/bin/awk

BEGIN {
	FS = "\t";
	file_idx = 0;
}

FNR == 1 {
	file_idx++;
}

# Seq and Struct Homology File
file_idx < 3 && FNR > 1 {
	vs_id_field = $2;
	n = split(vs_id_field, vs_id_array, ";");
	for(i = 1; i <= n; i++) {
		if(file_idx == 1) {
			seq[$1][vs_id_array[i]] = "TRUE";
			seq_score[$1][vs_id_array[i]] = $3;
		} else if(file_idx == 2) {
			struct[$1][vs_id_array[i]] = "TRUE";
			struct_score[$1][vs_id_array[i]] = $3;
		}
	}
}

file_idx == 3 && FNR == 1 {
	header = $0 \
		"\tSum of All Ratios" \
		"\tHomologous by Seq" \
		"\tSeq Alignment Score" \
		"\tHomologous by Struct" \
		"\tStruct Alignment Score"; 
	print header;
}

file_idx == 3 && FNR > 1 {
	prefix = $0;
	sum_score = $9 + $10 + $11 + $12 + $13 + $14 + $15 + $16;

	is_seq_hom = seq[$1][$3]; 
	cur_seq_score = seq_score[$1][$3]; 
	if(!is_seq_hom) { is_seq_hom = "FALSE"; }

	is_struct_hom = struct[$1][$3]; 
	cur_struct_score = struct_score[$1][$3]; 
	if(!is_struct_hom) { is_struct_hom = "FALSE"; }

	line = prefix "\t" sum_score \
		"\t" is_seq_hom "\t" cur_seq_score \
		"\t" is_struct_hom "\t" cur_struct_score;
	print line;
}