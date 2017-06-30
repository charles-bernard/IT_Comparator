#!/usr/bin/awk

function get_gene_name(attribute,  gene_name) {
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
		gene_name = "unknown n°" ug;
		ug++;
	}
	return gene_name;
}


function write_gene_tag(i) {
	tag = "gene_name=" gene_name[i] \
		";start=" gene_start[i] \
		";end=" gene_end[i];
	return tag;
}


function get_genomic_landscape(start, end, d) {
	# The env variable 'd' tells the function where 
	# to begin to read the gene dictionary

	# As long as the terminator is located after a gene end
	# Get the next gene.
	while(start > gene_end[d] && gene_end[d]) {
		d++;
	}

	# idx of first gene before IT
	bf = d - 1;
	bf_strand = gene_strand[bf];
	# Counter for nb of gene(s) before IT on same strand
	nb_bf = 1;
	# Get info on these series of genes
	while(gene_strand[bf] == bf_strand) {
		bf_series[nb_bf] = write_gene_tag(bf);
		nb_bf++;
		bf--;
	}
	# Concatenate these info into a tag
	bf_series_tag = bf_series[nb_bf-1];
	for(i = nb_bf-2; i >= 1; i--) {
		bf_series_tag = bf_series_tag "|" bf_series[i];
	}
	# Get gene tag of gene breaking the series
	bf_break = write_gene_tag(bf);
	delete bf_series;

	# Apply the same method for after series
	af = d;
	while(gene_start[af] < end) {
		af++;
	}
	af_strand = gene_strand[af];
	nb_af = 1;
	while(gene_strand[af] == af_strand) {
		af_series[nb_af] = write_gene_tag(af);
		nb_af++;
		af++;
	}
	af_series_tag = af_series[1];
	for(i = 2; i < nb_af; i++) {
		af_series_tag = af_series_tag "|" af_series[i];
	}
	af_break = write_gene_tag(af);
	delete af_series;

	genomic_landscape = bf_strand \
		"\t" nb_bf-1 \
		"\t" bf_series_tag \
		"\t" bf_break \
		"\t" af_strand \
		"\t" nb_af-1 \
		"\t" af_series_tag \
		"\t" af_break;

	return genomic_landscape;
}

BEGIN {
	FS = "\t";
	file_idx = 0;

	# gene counter
	g = 0;
	# counter for unknown gene
	ug = 0;
	# terminator counter
	t = 0;
	# counter for index of gene dictionary;
	# meant to avoid redundant reading;
	d = 0;
}

FNR == 1 {
	file_idx++;
}


# The input file 1 is the annotation of the genome
# Lines starting with a comment are ignored
file_idx == 1 && !/^#.*$/ {
	region = $3;
	# Only the Genes/Coding Sequences matter for us
	if(region == "gene") {
		# Important Fields
		gene_start[g] = $4;
		gene_end[g] = $5;
		gene_strand[g] = $7;
		gene_name[g] = get_gene_name($9);
		g++;
	}
}

# The input file 2 is the list of terminators
file_idx == 2 && FNR == 1 {
	if(NF == 6) {
		header = $0 \
			"\tStrand First Gene BEFORE IT" \
			"\tLength BEFORE series" \
			"\tSeries of Gene(s) BEFORE IT on Same Strand" \
			"\tGene on Opp. Strand Breaking BEFORE Series" \
			"\tStrand First Gene AFTER IT" \
			"\tLength AFTER series" \
			"\tSeries of Gene(s) AFTER IT on Same Strand" \
			"\tGene on Opp. Strand Breaking AFTER Series"; 
	} else {
		header = $0;
	}
}

# Only IT with no genomic landscape assigned will be annotated
file_idx == 2 && FNR > 1 && NF == 6 {
	prefix = $0;
	start = $1;
	end = $2;
	# The fct get_genomic_landscape will increment the dictionary index 'd'
	suffix = get_genomic_landscape(start, end, d);
	#print gene_start[0]
	line[t] = prefix "\t" suffix;
	t++;
}

file_idx == 3 && FNR > 1 && NF >= 12 {
	line[t] = $0; t++;
}

END {
	printf("%s\n", header);
	for(i = 0; i < t; i++) {
		printf("%s\n", line[i]);
	}
}