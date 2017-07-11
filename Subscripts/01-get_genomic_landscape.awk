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


function get_genomic_landscape(t, d,  i) {
	# The env variable 'd' tells the function where 
	# to begin to read the gene dictionary

	# As long as the terminator is located after a gene end
	# Get the next gene.
	while(it_start[t] > gene_end[d] && gene_end[d]) {
		d++;
	}

	# idx of first gene before IT
	g_bf = d - 1;
	g_bf_strand = gene_strand[g_bf];

	# previous IT on same strand
	it_bf = t-1; 
	while(it_strand[it_bf] && it_strand[it_bf] != g_bf_strand) { it_bf--; }
	if(!it_strand[it_bf]) { it_end[it_bf] = 0; }

	# Counter for nb of gene(s) before IT on same strand
	nb_g_bf = 1;
	# Get info on these series of genes
	while(gene_strand[g_bf] == bf_strand && gene_start[g_bf] > it_end[it_bf]) {
		bf_series[nb_g_bf] = write_gene_tag(g_bf);
		nb_g_bf++;
		g_bf--;
	}
	# Concatenate these info into a tag
	bf_series_tag = bf_series[nb_g_bf-1];
	for(i = nb_g_bf-2; i >= 1; i--) {
		bf_series_tag = bf_series_tag "|" bf_series[i];
	}

	# Get tag of feature breaking the series
	if(gene_start[g_bf] <= it_end[it_bf]) {
		bf_break = "cause=IT;start=" it_start[it_bf] ";end=" it_end[it_bf]; 
	} else {
		bf_break = "cause=Opposite_gene;" write_gene_tag(g_bf);
	}
	delete bf_series;

	# Apply the same method for after series
	g_af = d;
	g_af_strand = gene_strand[g_af];
	
	it_af = t+1; 
	while(it_strand[it_af] && it_strand[it_af] != af_strand) { it_af++; }
	if(!it_strand[it_af]) { it_start[it_af] = 10^9; }

	nb_g_af = 1;
	while(gene_strand[g_af] == af_strand && gene_end[g_af] < it_start[it_af]) {
		af_series[nb_g_af] = write_gene_tag(g_af);
		nb_g_af++;
		g_af++;
	}
	af_series_tag = af_series[1];
	for(i = 2; i < nb_g_af; i++) {
		af_series_tag = af_series_tag "|" af_series[i];
	}
	if(gene_end[g_af] >= it_start[it_af]) {
		af_break = "cause=IT;start=" it_start[it_af] ";end=" it_end[it_af]; 
	} else {
		af_break = "cause=Gene;" write_gene_tag(g_af);
	}
	delete af_series;

	genomic_landscape = bf_strand \
		"\t" nb_g_bf-1 \
		"\t" bf_series_tag \
		"\t" bf_break \
		"\t" af_strand \
		"\t" nb_g_af-1 \
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
			"\tFeature Breaking BEFORE Series" \
			"\tStrand First Gene AFTER IT" \
			"\tLength AFTER series" \
			"\tSeries of Gene(s) AFTER IT on Same Strand" \
			"\tFeature Breaking AFTER Series"; 
	} else {
		header = $0;
	}
}

# Only IT with no genomic landscape assigned will be annotated
file_idx == 2 && FNR > 1 && NF == 6 {
	prefix[t] = $0;
	it_start[t] = $1;
	it_end[t] = $2;
	it_strand[t] = $3;
	t++;
}

file_idx == 3 && FNR > 1 && NF > 6 {
	line[t] = $0; t++;
}

END {
	printf("%s\n", header);
	for(i = 0; i < t; i++) {
		# The fct get_genomic_landscape will increment 
		# the dictionary index 'd'
		suffix = get_genomic_landscape(i, d);
		line = prefix[i] "\t" suffix;
		printf("%s\n", line);
	}
}