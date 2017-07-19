#!/usr/bin/awk

BEGIN {
	FS = "\t";

	# counter for gene
	g = 0;

	# threshold
	th = 0.5;
}

NR > 1 {
	ref[g] = $1;
	score = $3;
	if(score < th) {
		hom[g] = $2;
		idx_hom[hom[g]] = g;
	} else {
		hom[g] = ".";
	}
	

	if(!occured_hom[hom[g]] && hom[g] != ".") {
		best_score[hom[g]] = score;
		best_idx_hom[hom[g]] = g; 
		occured_hom[hom[g]] = 1;
	} else if(occured_hom[hom[g]] && hom[g] != ".") {
		if(score < best_score[hom[g]]) {
			hom[best_idx_hom[$2]] = ".";
			best_score[hom[g]] = score;
			best_idx_hom[hom[g]] = g;
		} else {
			hom[g] = ".";
		}
	}

	g++;
}

END {
	printf("RefName\tHomName\n");
	for(i = 0; i < g; i++) {
		printf("%s\t%s\n", ref[i], hom[i]);
	}
}