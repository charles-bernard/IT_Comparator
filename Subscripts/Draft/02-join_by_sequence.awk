#!/usr/bin/awk

function abs(x) {
	return x < 0 ? -x : x;
}


function correct_and_normalize(alg_score) {
	# We don't want short IT matching longer IT get too much penalty
	# just because of a differential of length
	# gaps necessarily introduced because of a differential of length 
	# will have therefore only a cost of -1 instead of -2
	alg_score = alg_score + abs(ref_len - cur_len);

	# Normalize score
	if(ref_len > cur_len) {
		lg_len = ref_len;
	} else {
		lg_len = cur_len;
	}
	worst_case = lg_len * gap_cost;
	best_case = lg_len;
	alg_score = (alg_score - worst_case)/(best_case - worst_case);
	return alg_score;
}


function NW_matrix(ref_seq, cur_seq) {
	# This function builds the Needleman-Wunsch Score Matrix
	delete f;

	# Initialize matrix
	for(i = 0; i <= ref_len; i++) {
		f[i][0] = gap_cost * i;
	}
	for(j = 0; j <= cur_len; j++) {
		f[0][j] = gap_cost * j;
	}
	# Fill matrix
	for(i = 1; i <= ref_len; i++) {
		for(j = 1; j <= cur_len; j++) {
			ref_nt = substr(ref_seq, i, 1);
			comp_nt = substr(cur_seq, j, 1);

			match_score = f[i-1][j-1] + s[ref_nt][comp_nt];
			delete_score = f[i-1][j] + gap_cost; 
			insert_score = f[i][j-1] + gap_cost;

			if(match_score >= delete_score && match_score >= insert_score) {
				f[i][j] = match_score;
			} else if(delete_score > match_score && delete_score >= insert_score) {
				f[i][j] = delete_score;
			} else {
				f[i][j] = insert_score;
			}
		}
	}
}


function NW_alg(ref_seq, cur_seq,  \
	alg_score) {
	# Needleman-Wunsch alignment
	NW_matrix(ref_seq, cur_seq);

	i = ref_len; 
	j = cur_len;
	ref_alg = cur_alg = "";
	alg_score = 0;

	while(i > 0 || j > 0) {
		if(i > 0) { ref_nt = substr(ref_seq, i, 1); }
		if(j > 0) { cur_nt = substr(cur_seq, j, 1); }

		if(i > 0 && j > 0 && f[i][j] == f[i-1][j-1] + s[ref_nt][cur_nt]) {
			ref_alg = ref_nt ref_alg;
			cur_alg = cur_nt cur_alg;
			i--;
			j--;
			alg_score += s[ref_nt][cur_nt];
		} else if(i > 0 && f[i][j] == f[i-1][j] + gap_cost) {
			ref_alg = ref_nt ref_alg;
			cur_alg = "-" cur_alg;
			i--;
			alg_score += gap_cost;
		} else {
			ref_alg = "-" ref_alg;
			cur_alg = cur_nt cur_alg;
			j--;
			alg_score += gap_cost;
		}
	}
	alg_score = correct_and_normalize(alg_score);

	return alg_score ";" ref_alg "\\\\" cur_alg;
}


function best_match(ref_seq, seq, n,  i) {
	ref_len = length(ref_seq);
	for(i = 0; i < n; i++) {
		cur_seq = seq[i];
		cur_len = length(cur_seq)
		alg_out[i] = NW_alg(ref_seq, cur_seq);
		split(alg_out[i], alg_out_fields, ";");
		alg_score[i] = alg_out_fields[1];
		alg[i] = alg_out_fields[2];
	}
	
	asorti(alg_score, sorted_order, "@val_num_desc");

	best_idx = sorted_order[1];
	best_alg_score = alg_score[best_idx];
	best_alg = alg[best_idx];

	if(best_alg_score == alg_score[sorted_order[2]]) {
		print "#############"
		print "DUPLICATE!!!!"
		print "#############"
	}

	delete sorted_order;
	delete alg_out;
	delete alg_score;
	delete alg;

	return best_idx ";" best_alg_score ";" best_alg;
}


BEGIN {
	FS = "\t";
	file_idx = 0;
	# counter for nb of ITs to compare with reference IT
	c = 0;

	# by is sent to this script
	if(by == "sequence") {
		seq_field = 4;
	} else if(by == "structure") {
		seq_field = 5;
	}

	# Alignement Scores
	# For primary sequences
	s["A"]["A"] = s["T"]["T"] = s["C"]["C"] = s["G"]["G"] = 1;
	s["A"]["T"] = s["T"]["A"] = s["C"]["G"] = s["G"]["C"] = 0.25;
	s["A"]["C"] = s["A"]["G"] = s["T"]["C"] = s["T"]["G"] = -1;
	s["C"]["A"] = s["C"]["T"] = s["G"]["A"] = s["G"]["T"] = -1;
	# For secondary structures
	s["("]["("] = s[")"][")"] = s["."]["."] = 1;
	s["("]["."] = s["."]["("] = s[")"]["."] = s["."][")"] = -1;
	s["("][")"] = s[")"]["("] = -2;
	# cost for gap
	gap_cost = -2;
}


FNR == 1 {
	file_idx++;
}


# File to compare with Reference List
file_idx == 1 && FNR > 1 {
	seq[c] = $seq_field;
	line[c] = $0;
	c++;
}


file_idx == 2 && FNR == 1 {
	header = $0 "\tAlignement Score\tAlignement\t" $0;
	printf("%s\n", header);
}


file_idx == 2 && FNR > 1  {
	ref_seq = $seq_field;
	ref_line = $0;

	best_out = best_match(ref_seq, seq, c);
	split(best_out, best_out_fields, ";");
	best_idx = best_out_fields[1];
	best_alg_score = best_out_fields[2];
	best_alg = best_out_fields[3];
	#print FNR-1 "/" c ;
	
	printf("%s\t%s\t%s\t%s\n", ref_line, best_alg_score, best_alg, line[best_idx]);
}
