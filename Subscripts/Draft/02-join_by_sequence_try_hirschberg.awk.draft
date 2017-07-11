#!/usr/bin/awk

function abs(x) {
	return x < 0 ? -x : x;
}


function rev(seq,  \
	i,  nt,  rev_seq) {
	rev_seq = "";
	for(i = length(seq); i > 0; i--) {
		nt = substr(seq, i, 1);
		rev_seq = rev_seq nt;
	}
	return rev_seq;
}


function NW_matrix(ref_seq, cur_seq,  \
	ref_len,  cur_len,  i,  j,  ref_nt,  cur_nt,  \
	match_score,  delete_score,  insert_score) {
	# This function builds the Needleman-Wunsch Score Matrix 'f' (global variable)
	delete f;
	ref_len = length(ref_seq);
	cur_len = length(cur_seq);

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
			cur_nt = substr(cur_seq, j, 1);

			match_score = f[i-1][j-1] + s[ref_nt][cur_nt];
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


function NW_score(ref_seq, cur_seq,  \
	ref_len,  cur_len,  last_line,  j) {
	# returns the last line of the NW score matrix
	ref_len = length(ref_seq);
	cur_len = length(cur_seq);

	NW_matrix(ref_seq, cur_seq);

	last_line = f[ref_len][0];
	for(j = 1; j <= cur_len; j++) {
		last_line = last_line ";" f[ref_len][j];
	}

	return last_line;
}


function NW_alg(ref_seq, cur_seq,  \
	i,  j,  ref_len,  cur_len,  alg_score,   \
	ref_alg,  cur_alg,  ref_nt,  cur_nt) {
	# Needleman-Wunsch alignment
	NW_matrix(ref_seq, cur_seq);

	i = ref_len = length(ref_seq); 
	j = cur_len = length(cur_seq);
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

	return alg_score "\n" ref_alg "\n" cur_alg;
}


function min_NW_LR(NW_L, NW_R,  \
	L_n,  R_n,  scoreL,  scoreR,  i,  j,  max,  cur_mid) {
	L_n = split(NW_L, scoreL, ";");
	R_n = split(NW_R, scoreR, ";");
	if(L_n > R_n) { n = L_n; } else { n = R_n; }

	j = R_n;
	max = -1000;
	cur_mid = 0;

	for(i = 1; i <= n; i++) {
		if(scoreL[i] + scoreR[j] > max) {
			max = scoreL[i] + scoreR[j];
			cur_mid = i-1;
		}
		j--;
	}

	return cur_mid;
}


function Hirschberg_alg(ref_seq, cur_seq,  \
	ref_alg,  cur_alg,  ref_len,  cur_len,  i,  j,  \
	ref_mid,  cur_mid,  NW_L,  NW_R,  HL,  RL) {
	ref_alg = cur_alg = "";
	ref_len = length(ref_seq);
	cur_len = length(cur_seq);

	if(ref_len == 0) {
		for(j = 1; j <= cur_len; j++) {
			ref_alg = ref_alg "-";
			cur_alg = cur_alg substr(cur_seq, j, 1);
		}
		alg_score = cur_len * gap_cost;
	} else if(cur_len == 0) {
		for(i = 1; i <= ref_len; i++) {
			ref_alg = ref_alg substr(ref_seq, i, 1);
			cur_alg = cur_alg "-";
		}
		alg_score = ref_len * gap_cost;
	} else if(ref_len == 1 || cur_len == 1) {
		NW_out = NW_alg(ref_seq, cur_seq);
		split(NW_out, NW_fields, "\n");
		alg_score = NW_fields[1];
		ref_alg = NW_fields[2];
		cur_alg = NW_fields[3];
		delete NW_fields;
	} else {
		ref_mid = int(ref_len/2);

		NW_L = NW_score(substr(ref_seq, 1, ref_mid), cur_seq);
		NW_R = NW_score(rev(substr(ref_seq, ref_mid+1, ref_len-ref_mid)), rev(cur_seq));
		cur_mid = min_NW_LR(NW_L, NW_R);

		Hirschberg_alg(substr(ref_seq, 1, ref_mid), 
			substr(cur_seq, 1, cur_mid));
		Hirschberg_alg(substr(ref_seq, ref_mid+1, ref_len-ref_mid), 
			substr(cur_seq, cur_mid+1, cur_len-cur_mid));
	}
	ref_after_alg = ref_after_alg ref_alg;
	cur_after_alg = cur_after_alg cur_alg;
	final_alg_score += alg_score;
}


function correct_and_normalize(ref_seq, cur_seq, alg_score) {
	ref_len = length(ref_seq);
	cur_len = length(cur_seq);
	# We don't want short IT matching longer IT get too much penalty
	# just because of a differential of length
	# gaps necessarily introduced because of a differential of length 
	# will have therefore only a cost of -1 instead of -2
	alg_score += abs(ref_len - cur_len);

	# Normalize score
	if(ref_len > cur_len) {
		lg_len = ref_len;
	} else {
		lg_len = cur_len;
	}
	worst_case = lg_len * gap_cost;
	best_case = lg_len * 2;
	print "alg_score " alg_score
	print "best_case " best_case
	print "worst_case" worst_case
	alg_score = (alg_score - worst_case)/(best_case - worst_case);
	print "alg_score " alg_score
	if (alg_score > 1) { print "difheurfhuyf"}
	return alg_score;
}


function best_match(ref_seq, seq, n,  \
	i,  cur_seq,  alg_score,  alg) {
	for(i = 0; i < n; i++) {
		cur_seq = seq[i];

		ref_after_alg = cur_after_alg = "";
		final_alg_score = 0;

		# the function will update the three global
		# variables above
		Hirschberg_alg(ref_seq, cur_seq);

		# print "-------------"
		# print ref_after_alg
		# print cur_after_alg
		# print "score=" final_alg_score
		# print "-------------"

		alg_score[i] = correct_and_normalize(ref_seq, cur_seq, final_alg_score);
		alg[i] = ref_after_alg "\\" cur_after_alg;

		# print "-------------"
		# print alg[i]
		# print "norm_score=" alg_score[i]
		# print "-------------"
	}
	
	asorti(alg_score, sorted_order, "@val_num_desc");

	best_idx = sorted_order[1];
	best_alg_score = alg_score[best_idx];
	best_alg = alg[best_idx];

	delete sorted_order;
	delete alg_score;
	delete alg;

	return best_idx "\t" best_alg_score "\t" best_alg;
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
	s["A"]["A"] = s["T"]["T"] = s["C"]["C"] = s["G"]["G"] = 2;
	s["A"]["T"] = s["T"]["A"] = s["C"]["G"] = s["G"]["C"] = -1;
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
	split(best_out, best_out_fields, "\t");
	best_idx = best_out_fields[1];
	best_alg_score = best_out_fields[2];
	best_alg = best_out_fields[3];
	print FNR "/" c ;
	
	printf("%s\t%s\t\t%s%s\n", ref_line, best_alg_score, best_alg, line[best_idx]);
}
