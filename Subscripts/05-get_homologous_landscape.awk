#!/usr/bin/awk

# Sorry the code is a bit messy, I coded that in a hurry...
# I will try to rearrange it soon ;).

function get_gene_name(attribute) {
	match(attribute, /^gene_name[= ]+"*[a-zA-Z0-9?.+:_-]+[ "]*;/, gene_tag);
	gsub(/gene_name[= ]+"*/, "", gene_tag[0]);
	gsub(/[ "]*;$/, "", gene_tag[0]);
	gene_name = gene_tag[0];
	delete gene_tag;
	return gene_name;
}

function get_break_event(attribute) {
	match(attribute, /^cause[= ]+"*[a-zA-Z0-9?.+:_-]+[ "]*;/, break_tag);
	gsub(/cause[= ]+"*/, "", break_tag[0]);
	gsub(/[ "]*;$/, "", break_tag[0]);
	break_event = break_tag[0];
	delete break_tag;
	return break_event;
}

function series_to_gene_array(max_g, n, series, type, t) {
	if(n > max_g) { limit = max_g; } else { limit = n; }
	split(series, gene_attribute, "|");

	for(i = 1; i <= limit; i++) {
		if(type == "bf") {
			if(file_idx == 1) {
				bf_genes[t][i-1] = get_gene_name(gene_attribute[n+1-i]);
			} else {
				ref_bf_genes[i-1] = get_gene_name(gene_attribute[n+1-i]);
			}
		} else if(type == "af") {
			if(file_idx == 1) {
				af_genes[t][i-1] = get_gene_name(gene_attribute[i]);
			} else {
				ref_af_genes[i-1] = get_gene_name(gene_attribute[i]);
			}
		}
	}
	delete gene_attribute;
}

function convert_to_contrary_series(string,  char) {
	len = length(string);
	contrary = "";
	for(i = 1; i <= len; i++) {
		if(i == 1) {
			char = substr(string, len, 1);
		} else if(i == len) {
			char = substr(string, 1, 1);
		} else {
			char = substr(string, i, 1);
		}
		contrary = contrary char;
	}
	return contrary;
}

function build_reference_code(max_g, n, type, direction, break_event) {
	if(n > max_g) {
		break_event = max_g;
	}
	break_char = code[type][break_event];

	i = 0;
	code_string = "";
	while(i < n && i < max_g) {
		if(type == "bf") {
			code_string = code[type][i] code_string;
		} else {
			code_string = code_string code[type][i];
		}
		i++;
	}

	if(direction == "opposite") {
		code_string = tolower(code_string);
	} else if(break_event == "Opposite_gene") {
		break_char = tolower(break_char);
	}

	if(type == "bf") {
		code_string = break_char code_string code[type]["mainIT"];
	} else {
		code_string = code[type]["mainIT"] code_string break_char;
	}

	return code_string;
}

function initialize_vs_code(max_g, n, type, direction, break_event, t) {
	if(n > max_g) {
		break_event = max_g;
	}
	break_char = code[type][break_event];

	i = 0;
	code_string = cur_char = "";
	while(i < n && i < max_g) {
		if(type == "bf") {
			if(bf_genes[t][i] ~ /^H:.*/) {
				cur_char = mismatch_homolog;
			} else {
				cur_char = mismatch_unknown;
			}
			code_string = cur_char code_string;
		} else {
			if(af_genes[t][i] ~ /^H:.*/) {
				cur_char = mismatch_homolog;
			} else {
				cur_char = mismatch_unknown;
			}
			code_string = code_string cur_char;
		}
		i++;
	}

	if(direction == "opposite") {
		code_string = tolower(code_string);
	} else if(break_event == "Opposite_gene") {
		break_char = tolower(break_char);
	}

	if(type == "bf") {
		code_string = break_char code_string code[type]["mainIT"];
	} else {
		code_string = code[type]["mainIT"] code_string break_char;
	}

	return code_string;
}

function find_match(max_g, contrary_series, type, 
	ref_code, ref_id, ref_genes, ref_n, 
	vs_code, vs_id, vs_genes, vs_n, t) {
	id_list = score_list = code_list = present_list = \
	order_list = orientation_list = flipped_list = "";
	k = 0;

	if(ref_n > max_g) { ref_n = max_g; }

	for(cur_t = 0; cur_t < t; cur_t++) {

		nb_g_present = 0;
		nb_g_ordered = 0;
		nb_g_oriented = 0;

		cur_code_string = vs_code[cur_t];

		if(vs_n[cur_t] > max_g) { vs_n[cur_t] = max_g; }

		for(i = 0; i < ref_n; i++) { 
			ref_char = substr(ref_code, i+2, 1);
			cur_char = ""; 

			if(type == "af") {
				ref_idx = i;
			} else if(type == "bf" ) {
				ref_idx = ref_n - i - 1;
			}
			for(j = 0; j < vs_n[cur_t]; j++) {

				if(type == "af") {
					idx = j;
				} else if(type == "bf" ) {
					idx = vs_n[cur_t] - j - 1;
				}	

				if(("H:" ref_genes[ref_idx]) == vs_genes[cur_t][idx]) {

					nb_g_present++;
					
					if(i == j) {
						nb_g_ordered++;
					}
					
					cur_char = substr(cur_code_string, j+2, 1);
					if(cur_char ~ /[A-Z]/) {
						cur_char = code[type][ref_idx];
						if(ref_char ~ /[A-Z]/) {
							nb_g_oriented++;
						} 
					} else {
						cur_char = tolower(code[type][ref_idx]);
						if(ref_char ~ /[a-z]/) {
							nb_g_oriented++;
						} 
					}

					cur_code_string = substr(cur_code_string, 1, j+1) \
						cur_char substr(cur_code_string, j+3);
				}
			}
		}

		if(nb_g_present) {
			
			if(contrary_series) {
				cur_code_string = convert_to_contrary_series(cur_code_string);
				cur_is_flipped = "TRUE";
			} else {
				cur_is_flipped = "FALSE";
			}

			if(k == 0) {
				id_list = id[cur_t];
				score_list = "(P=" nb_g_present "/" ref_n \
					",O=" nb_g_ordered "/" ref_n \
					",D=" nb_g_oriented "/" ref_n ")";;
				code_list = cur_code_string;
				present_list = nb_g_present / ref_n;
				order_list = nb_g_ordered / ref_n;
				orientation_list = nb_g_oriented / ref_n;
				flipped_list = cur_is_flipped;
			} else {
				id_list = id_list ";" id[cur_t];
				score_list = score_list ";(P=" nb_g_present "/" ref_n \
					",O=" nb_g_ordered "/" ref_n \
					",D=" nb_g_oriented "/" ref_n ")";;
				code_list = code_list ";" cur_code_string;
				present_list = present_list ";" nb_g_present / ref_n;
				order_list = order_list ";" nb_g_ordered / ref_n;
				orientation_list = orientation_list ";" nb_g_oriented / ref_n;
				flipped_list = flipped_list ";" cur_is_flipped;
			}
			k++;
		}
		
	}
	if(id_list) {
		return id_list "|" flipped_list "|" code_list "|" score_list \
			"|" present_list "|" order_list "|" orientation_list;
	}
	return "";
}

function assemble_matches(ref_code, bf_bf, bf_af, af_af, af_bf,
	  i,  j,  k) {
	matches_array[0] = bf_bf; matches_array[1] = bf_af;
	matches_array[2] = af_af; matches_array[3] = af_bf;

	i_bf_bf = 0; i_bf_af = 1; i_af_af = 2; i_af_bf = 3;

	j_id = 0; j_flip = 1; j_code = 2; j_score = 3;
	j_present = 4; j_order = 5; j_orientation = 6;

	# get attributes of each IT
	for(i = 0; i < 4; i++) {
		natts = split(matches_array[i], matches_attributes, "|");
		for(j = 0; j < natts; j++) {
			ntags[i] = split(matches_attributes[j+1], tmp_tags, ";");
			for(k = 0; k < ntags[i]; k++) {
				tag[i][j][k] = tmp_tags[k+1];
			}
			delete tmp_tags;
		}
		delete matches_attributes;
	}
	delete matches_array;

	# see if context is preserved on both side of the IT
	for(kA = 0; kA < ntags[i_bf_bf]; kA++) {
		cur_a_id = tag[i_bf_bf][j_id][kA];
		k_other_side[cur_a_id] = "";
		for(kB = 0; kB < ntags[i_af_af]; kB++) {
			cur_b_id = tag[i_af_af][j_id][kB]
			if(cur_a_id == cur_b_id) {
				k_other_side[cur_a_id] = kB;
			}
		}
	}
	for(kA = 0; kA < ntags[i_bf_af]; kA++) {
		cur_a_id = tag[i_bf_af][j_id][kA];
		k_other_side[cur_a_id] = "";
		for(kB = 0; kB < ntags[i_af_bf]; kB++) {
			cur_b_id = tag[i_af_bf][j_id][kB]
			if(cur_a_id == cur_b_id) {
				k_other_side[cur_a_id] = kB;
			}
		}
	}

	# Build final tags
	line = "";
	for(i = 0; i < 4; i++) {
		for(k = 0; k < ntags[i]; k++) {
			cur_id = tag[i][j_id][k];
			cur_line = id_ = flip = "";
			bf_code_ = bf_score = bf_present = bf_order = bf_orientation = "";
			af_code_ = af_score = af_present = af_order = af_orientation = "";
			if(!already_visited[cur_id]) {
				id_ = cur_id;
				flip = tag[i][j_flip][k];
				if(i == 0 || i == 1) {
					bf_code_ = tag[i][j_code][k];
					bf_score = tag[i][j_score][k];
					bf_present = tag[i][j_present][k];
					bf_order = tag[i][j_order][k];
					bf_orientation = tag[i][j_orientation][k];
					if(k_other_side[cur_id] != "") {
						if(i == 0) { i_other_side = 2; } else { i_other_side = 3; }
						af_code_ = tag[i_other_side][j_code][k_other_side[cur_id]];
						af_score = tag[i_other_side][j_score][k_other_side[cur_id]];
						af_present = tag[i_other_side][j_present][k_other_side[cur_id]];
						af_order = tag[i_other_side][j_order][k_other_side[cur_id]];
						af_orientation = tag[i_other_side][j_orientation][k_other_side[cur_id]];
					}
				} else if(i == 2 || i == 3) {
					af_code_ = tag[i][j_code][k];
					af_score = tag[i][j_score][k];
					af_present = tag[i][j_present][k];
					af_order = tag[i][j_order][k];
					af_orientation = tag[i][j_orientation][k];
				}
			}
			already_visited[cur_id] = 1;

			cur_line = ref_id "\t" ref_bf_code ref_af_code "\t" \
				id_ "\t" flip "\t" bf_code_ "\t" af_code_ "\t" bf_score "\t" af_score "\t" \
				bf_present "\t" bf_order "\t" bf_orientation "\t" \
				af_present "\t" af_order "\t" af_orientation;
			if(!line && id_) {
				line = cur_line;
			} else if(id_) {
				line = line "\n" cur_line;
			}
		}
	}
	delete k_other_side; 
	delete tag;
	delete ntags;
	delete already_visited;
	if(line) {
		return line;
	}
	return ref_id "\t" ref_bf_code ref_af_code "\t\t\t\t\t\t\t\t\t\t\t\t"
}


BEGIN {
	FS = "\t";
	file_idx = 0;

	# counter for IT
	t = 0;

	header = "Ref. ID" \
		"\tRef. Code" \
		"\tVs ID" \
		"\tis Vs Flipped?" \
		"\tVs Bf Code" \
		"\tVs Af Code" \
		"\tVs Bf Score" \
		"\tVs Af Score" \
		"\tVs Bf NgPresent/totNg" \
		"\tVs Bf NgOrdered/totNg" \
		"\tVs Bf NgOriented/totNg" \
		"\tVs Af NgPresent/totNg" \
		"\tVs Af NgOrdered/totNg" \
		"\tVs Af NgOriented/totNg";
	printf("%s\n", header);

	max_g = 5;

	# Code for landscape context:
	code["bf"][0] = "A";
	code["bf"][1] = "B";
	code["bf"][2] = "C";
	code["bf"][3] = "D";
	code["bf"][4] = "E";

	code["af"][0] = "Z";
	code["af"][1] = "Y";
	code["af"][2] = "X";
	code["af"][3] = "W";
	code["af"][4] = "V";

	code["bf"]["Opposite_gene"] = "K";
	code["af"]["Opposite_gene"] = "Q";

	code["bf"][max_g] = code["af"][max_g] = ".";

	code["bf"]["IT"] = code["af"]["IT"] = ":";
	code["bf"]["mainIT"] = code["af"]["mainIT"] = "!";

	mismatch_homolog = "H";
	mismatch_unknown = "N";

	direction["+"]["bf"] = direction["-"]["af"] = "towards";
	direction["+"]["af"] = direction["-"]["bf"] = "opposite";

	# Rules:
	# if gene towards
	#	-> toupper(Code_Letter);
	# else:
	#	-> tolower(Code_Letter);

	# Example:
	# Reference Code:
	# :BA!vwxQ
	# (IT!)-[ B }-[ A }-(IT!)-[ v }-[ w }-[ x }-{ Q ]
	# Vs code:
	# kHB!vwhX
	# { k ]-[ H }-[ B }-(IT!)-[ v }-[ w }-[ h }-{ X ]

}

FNR == 1 {
	file_idx++;
}

file_idx == 1 && FNR > 1 {
	id[t] = $1;

	bf_direction[t] = direction[$8]["bf"];
	n_g_bf[t] = $9;
	bf_series_tag = $10;
	bf_break_event[t] = get_break_event($11);

	# The function below fills the array bf_genes[t][@]
	series_to_gene_array(max_g, n_g_bf[t], bf_series_tag, "bf", t);

	bf_code[t] = initialize_vs_code(max_g, n_g_bf[t], 
		"bf", bf_direction[t], bf_break_event[t], t);

	af_direction[t] = direction[$12]["af"];
	n_g_af[t] = $13;
	af_series_tag = $14;
	af_break_event[t] = get_break_event($15);

	# The function below fill the array af_genes[t][@]
	series_to_gene_array(max_g, n_g_af[t], af_series_tag, "af", t);

	af_code[t] = initialize_vs_code(max_g, n_g_af[t], 
		"af", af_direction[t], af_break_event[t], t);

	t++;
}

file_idx == 2 && FNR > 1 {
	ref_id = $1;

	ref_bf_direction = direction[$8]["bf"];
	ref_n_g_bf = $9;
	ref_bf_series_tag = $10;
	ref_bf_break_event = get_break_event($11);

	# The function below fill the array ref_bf_genes[t][@]
	series_to_gene_array(max_g, ref_n_g_bf, ref_bf_series_tag, "bf", 0);

	ref_bf_code = build_reference_code(max_g, ref_n_g_bf, 
		"bf", ref_bf_direction, ref_bf_break_event);

	ref_af_direction = direction[$12]["af"];
	ref_n_g_af = $13;
	ref_af_series_tag = $14;
	ref_af_break_event = get_break_event($15);

	# The function below fill the array ref_af_genes[t][@]
	series_to_gene_array(max_g, ref_n_g_af, ref_af_series_tag, "af", 0);

	ref_af_code = build_reference_code(max_g, ref_n_g_af, 
		"af", ref_af_direction, ref_af_break_event);

	## FIND MATCHES ##
	match_bf_bf = find_match(max_g, 0, "bf", 
		ref_bf_code, ref_id, ref_bf_genes, ref_n_g_bf, 
		bf_code, id, bf_genes, n_g_bf, t);

	match_bf_af = find_match(max_g, 1, "bf", 
		ref_bf_code, ref_id, ref_bf_genes, ref_n_g_bf, 
		af_code, id, af_genes, n_g_af, t);

	match_af_af = find_match(max_g, 0, "af",
		ref_af_code, ref_id, ref_af_genes, ref_n_g_af, 
		af_code, id, af_genes, n_g_af, t);

	match_af_bf = find_match(max_g, 1, "af", 
		ref_af_code, ref_id, ref_af_genes, ref_n_g_af, 
		bf_code, id, bf_genes, n_g_bf, t);


	line = assemble_matches(ref_bf_code ref_af_code,
		match_bf_bf, match_bf_af,
		match_af_af, match_af_bf);

	print line;

	delete ref_bf_genes; delete ref_af_genes;
	ref_bf_genes[0] = ref_af_genes[0] = "";
}