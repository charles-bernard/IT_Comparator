#!/usr/bin/awk

function get_gene_name(attribute) {
	# The Gene tag is the match of the Regex 
	# against the 'attribute' string
	match(attribute, /; *(gene_)?[nN]ame[= ]+"*[a-zA-Z0-9.+:_-]+[ "]*;/, gene_tag);
	if(!gene_tag[0]) {
		# If gene has no name, get its ID
		match(attribute, /; *(gene_)?[iI]d[= ]+[a-zA-Z0-9.+:_-]+[ "]*;/, gene_tag);
	}
	gsub(/^; *(gene_)?[nNiI](ame|d)[= ]+"*/, "", gene_tag[0]);
	gsub(/[ "]*;$/, "", gene_tag[0]);
	gene_name = gene_tag[0];
	if(!gene_name) {
		gene_name = "unknown_n°" ug;
		ug++;
	}
	return gene_name;
}

function is_rna(attribute) {
	match(attribute, /^.*(CDS|protein_coding).*$/, cds);
	if(cds[0]) {
		return 0;
	}
	match(attribute, /^.*(RNA|rna).*$/, rna);
	if(!rna[0]) {
		return 0;
	} else {
		return 1;
	}
}

function write_seq(strand, start, end,  seq,  rev_seq) {
	# The terminator sequence is just a substring of the genome sequence
	seq = substr(genome_seq, start, end-start+1);

	if(strand == "+") {
		return seq;
	} else {
		# If the terminator is on the "-" strand
		# the sequence needs to be reversed
		rev_seq = "";

		# The sequence is read from the end to the start
		for(i = length(seq); i > 0; i--) {
			# Then each nt will be reversed
			# and concatenated to the previous
			nt = substr(seq, i, 1);

			# Check if nt is A,C,T or G 
			if(nt ~ /[ATCG]/) {
				rev_seq = rev_seq rev[nt];
			# A gap character "-" will be use if not
			} else {
				rev_seq = rev_seq "-";
			}
		}
		return rev_seq;
	}
}


function translate(seq,  protein) {
	n = length(seq);
	protein = "";
	nbAbnormalities = 0;

	for(i = 1; i < n; i += 3) {
		cur_aa = aa[substr(seq, i, 3)];

		if(cur_aa) {
			protein = protein cur_aa;
		} else {
			protein = protein "-";
			nbAbnormalities++;
		}

		if(i == 1 && cur_aa != "M") {
			# if first aa is not a methionine
			nbAbnormalities++;
		} else if(i + 3 < n && cur_aa == "*") {
			# if a stop codon occurs before the end
			nbAbnormalities++;
		} else if(i + 3 >= n && cur_aa != "*") {
			# if the seq is not terminated by a stop codon
			nbAbnormalities++;
		}
		
	}

	if(nbAbnormalities > 2) {
		return "abnormal_CDS\t" protein;
	} else {
		return "CDS\t" protein;
	}
}

BEGIN {
	FS = "\t";

	file_idx = 0;
	
	genome_seq = "";

	# Genetic Code
	aa["TAA"] = aa["TAG"] = aa["TGA"] = "*";
	aa["GCT"] = aa["GCC"] = aa["GCA"] = aa["GCG"] = "A";
	aa["TGT"] = aa["TGC"] = "C";
	aa["GAT"] = aa["GAC"] = "D";
	aa["GAA"] = aa["GAG"] = "E";
	aa["TTT"] = aa["TTC"] = "F";
	aa["GGT"] = aa["GGC"] = aa["GGA"] = aa["GGG"] = "G";
	aa["CAT"] = aa["CAC"] = "H";
	aa["ATT"] = aa["ATC"] = aa["ATA"] = "I";
	aa["AAA"] = aa["AAG"] = "K";
	aa["TTA"] = aa["TTG"] = aa["CTT"] = aa["CTC"] = aa["CTA"] = aa["CTG"] = "L"; 
	aa["ATG"] = "M";
	aa["AAT"] = aa["AAC"] = "N";
	aa["CCT"] = aa["CCC"] = aa["CCA"] = aa["CCG"] = "P";
	aa["CAA"] = aa["CAG"] = "Q";
	aa["CGT"] = aa["CGC"] = aa["CGA"] = aa["CGG"] = aa["AGA"] = aa["AGG"] = "R";
	aa["TCT"] = aa["TCC"] = aa["TCA"] = aa["TCG"] = aa["AGT"] = aa["AGC"] = "S";
	aa["ACT"] = aa["ACC"] = aa["ACA"] = aa["ACG"] = "T";
	aa["GTT"] = aa["GTC"] = aa["GTA"] = aa["GTG"] = "V";
	aa["TGG"] = "W";
	aa["TAT"] = aa["TAC"] = "Y";
	
	# NT complementation
	rev["A"] = "T";
	rev["T"] = "A";
	rev["C"] = "G";
	rev["G"] = "C";

	 # gene counter
	g = ug = 0;

	# Three files are sent as args to this script
	# - cds_file
	# - abn_cds_file
	# - nc_file

	header = "start" \
		"\tend" \
		"\tstrand" \
		"\tname" \
		"\tbiotype" \
		"\tsequence";
	printf("%s\n", header);
}

FNR == 1 {
	file_idx++;
}

# The input file 1 corresponds to the genome in fasta
file_idx == 1 && FNR > 1 {
	# Concatenate each line in the fasta to get
	# the whole genome sequence as a single string
	genome_seq = genome_seq $0;
}

# The input file 2 is the annotation of the genome
# Lines starting with a comment are ignored
file_idx == 2 && !/^#.*$/ {
	region = $3;
	# Only the Genes/Coding Sequences matter for us
	if(region == "gene") {
		# Important Fields
		start = $4;
		end = $5;
		strand = $7;
		seq = write_seq(strand, start, end);
		name = get_gene_name($9);
		isRNA = is_rna($9);
		if(isRNA) {
			seq_fields = "non_coding\t" seq;
		} else {
			seq_fields = translate(seq);
		}

		split(seq_fields, f, "\t");
		if(length(f[2]) < 4096) {
			printf("%s\t%s\t%s\t%s\t%s\n", 
				start, end, strand,
				name, seq_fields)
			if(f[1] == "CDS") {
				printf("%s\t%s\n", name, f[2]) > cds_file
			} else if(f[1] == "abnormal_CDS") {
				printf("%s\t%s\n", name, f[2]) > abn_cds_file
			} else {
				printf("%s\t%s\n", name, f[2]) > nc_file
			}
		} 
	}
}
