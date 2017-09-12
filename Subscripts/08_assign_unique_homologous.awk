#!/usr/bin/awk

function sort_dict_entry(dict, id, n,  
	i,  j,  cur_id) {
	# Simple insertion sort

	if(id ~ /IT[+][0-9]+/) { ref_strand = "+"; } else { ref_strand = "-"; }

	i = 2;

	while(i <= n) {

		cur_id = dict[id][i]["id"];
		if(cur_id ~ /IT[+][0-9]+/) { cur_strand = "+"; } else { cur_strand = "-"; }
		cur_c_score = dict[id][i]["c_score"];
		cur_match = dict[id][i]["match"];
		cur_a_score = dict[id][i]["a_score"];

		j = i - 1;
		better = 0;

		while(j > 0 && !better) {

			if(cur_match != "TRUE" && dict[id][j]["match"] == "TRUE") {
				better = 1;
			} else if(cur_match == dict[id][j]["match"]) {
				if(cur_c_score < dict[id][j]["c_score"]) {
					better = 1;
				} else if(cur_c_score == dict[id][j]["c_score"]) {
					if(dict[id][j]["id"] ~ /IT[+][0-9]+/) { comp_strand = "+"; } else { comp_strand = "-"; }
					if(ref_strand == comp_strand && ref_strand != cur_strand) {
						better = 1;
					} else if(cur_strand == comp_strand) {
						if(cur_a_score > dict[id][j]["a_score"]) {
							better = 1;
						}
					}
				}
			}
			if(!better) {
				dict[id][j+1]["id"] = dict[id][j]["id"];
				dict[id][j+1]["c_score"] = dict[id][j]["c_score"];
				dict[id][j+1]["match"] = dict[id][j]["match"];
				dict[id][j+1]["a_score"] = dict[id][j]["a_score"];
				j--;
			} 
		}

		dict[id][j+1]["id"] = cur_id;
		dict[id][j+1]["c_score"] = cur_c_score;
		dict[id][j+1]["match"] = cur_match;
		dict[id][j+1]["a_score"] = cur_a_score;
		i++;

	}
}

BEGIN {
	FS = "\t";

	t = 0;
	ref_id = "";
}

NR == 1 {
	header = $0;
	print header;
}

NR > 1 {
	if($1 == ref_id) {

		n_has[ref_id]++;
		vs_id = $3;
		if(n_is[vs_id]) {
			n_is[vs_id]++;
		} else {
			n_is[vs_id] = 1;
		}

	} else {
		
		if(ref_id) {
			sort_dict_entry(has_hom, ref_id, n_has[ref_id]);
		}
		list_it[t] = ref_id = $1;
		list_code[t] = $2;
		n_has[ref_id] = 0;
		if($3) { 
			vs_id = $3;
			n_has[ref_id] = 1;
			if(n_is[vs_id]) {
				n_is[vs_id]++;
			} else {
				n_is[vs_id] = 1;
			}
		} else {
			vs_id = "";
			n_has[ref_id] = 0;
		} 
		t++;
	}

	if(vs_id) {
		has_hom[ref_id][n_has[ref_id]]["id"] = $3;
		is_hom[vs_id][n_is[vs_id]]["id"] = $1;

		is_hom[vs_id][n_is[vs_id]]["c_score"] = \
			has_hom[ref_id][n_has[ref_id]]["c_score"] = $17;
		is_hom[vs_id][n_is[vs_id]]["match"] = \
			has_hom[ref_id][n_has[ref_id]]["match"] = $18;
		is_hom[vs_id][n_is[vs_id]]["a_score"] = \
			has_hom[ref_id][n_has[ref_id]]["a_score"] = $19;
		fields[ref_id][vs_id] = $0;
	} 


}

END {
	for(i = 0; i < t; i++) {
		cur_id = list_it[i];
		k = 1;
		assigned = 0;
		line = cur_id "\t" list_code[i] "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
		while(!assigned && k <= n_has[cur_id]) {
			vs_id = has_hom[cur_id][k]["id"];
			if(!has_been_sorted[vs_id]) {
				sort_dict_entry(is_hom, vs_id, n_is[vs_id]);
				has_been_sorted[vs_id] = 1;
			}
			if(cur_id == is_hom[vs_id][1]["id"]) {
				line = fields[cur_id][vs_id];
				assigned = 1;
			} else {
				k++;
			}
		}
		print line;
	}
}