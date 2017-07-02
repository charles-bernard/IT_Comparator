#!/usr/bin/awk

function alignment(ref_seq, compl_seq,  len_score,  i,  j) {
	ref_len = length(ref_seq);
	compl_len = length(compl_seq);

	if(ref_len >= compl_len) {
		lg_len = ref_len;
		sh_len = compl_len;
		lg_seq = ref_seq;
		sh_seq = compl_seq; 
	} else {
		lg_len = compl_len;
		sh_len = ref_len;
		lg_seq = compl_seq;
		sh_seq = ref_seq;
	}

	diff_len = lg_len - sh_len;
	len_score = 1/(diff_len + 1);

	best_alg_score = 0;
	for(i = 1; i <= diff_len; i++) {
		cur_alg_score = 0;
		for(j = 1; j < 1 + sh_len; j++) {
			ln_nt = substr(lg_seq, j, 1);
			sh_nt = substr(sh_seq, j-i+1, 1);
			cur_alg_score += s[ln_nt][sh_nt];
		}
		if(cur_alg_score > best_alg_score) {
			best_alg_score = cur_alg_score;
		}
	}
	return best_alg_score / sh_len  "\t" len_score;
}


function best_match(ref_seq, seq, n,  i) {
	ref_len = length(ref_seq);
	for(i = 0; i < n; i++) {
		cur_seq = seq[i];
		alg_out[i] = alignment(ref_seq, cur_seq);
		split(alg_out[i], f, "\t");
		alg_score[i] = f[1];
		len_score[i] = f[2];
	}
	
	asorti(alg_score, sorted_order, "@val_num_desc");
	best_alg = alg_score[sorted_order[1]];
	best_len = len_score[sorted_order[1]];
	best_idx = sorted_order[1];
	k = 2;
	while(alg_score[sorted_order[k]] == best_alg) {
		if(len_score[sorted_order[k]] > best_len) {
			best_idx = sorted_order[k];
			best_len = len_score[best_idx];
		}
		k++;
	}

	delete alg_out;
	delete alg_score; 
	delete len_score;
	delete sorted_order;

	return best_idx "\t" best_alg;
}


BEGIN {
	FS = "\t";
	file_idx = 0;
	# counter for nb of ITs to compare with reference IT
	c = 0;

	# Alignement Scores
	s["A"]["A"] = s["T"]["T"] = s["C"]["C"] = s["G"]["G"] = 1;
	s["A"]["T"] = s["T"]["A"] = s["C"]["G"] = s["G"]["C"] = 0.75;
	s["A"]["C"] = s["A"]["G"] = s["T"]["C"] = s["T"]["G"] = 0;
	s["C"]["A"] = s["C"]["T"] = s["G"]["A"] = s["G"]["T"] = 0;
}


FNR == 1 {
	file_idx++;
}


# File to compare with Reference List
file_idx == 1 && FNR > 1 {
	seq[c] = $4;
	line[c] = $0;
	c++;
}


file_idx == 2 && FNR == 1 {
	header = $0 "\tAlignement Score\t" $0;
	printf("%s\n", header);
}


file_idx == 2 && FNR > 1  {
	ref_seq = $4;
	ref_line = $0;

	best_out = best_match(ref_seq, seq, c);
	split(best_out, f, "\t");
	best_i = f[1];
	best_alg_score = f[2];
	
	printf("%s\t%s\t%s\n", ref_line, best_alg_score, line[best_i]);
}
