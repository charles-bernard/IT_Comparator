#!/usr/bin/awk

function write_gene_tag(i) {
	tag = "gene_name=" gene_name[i] \
		";start=" gene_start[i] \
		";end=" gene_end[i];
	return tag;
}


function get_genomic_landscape(t, d, breakIT,  i) {
	# The env variable 'd' tells the function where 
	# to begin to read the gene dictionary

	# As long as the terminator is located after a gene end
	# Get the next gene.
	while(it_start[t] > gene_end[d] && gene_end[d]) {
		d++;
	}

	## BEFORE SERIES ##

	# The index & the strand of the first gene before the IT
	g_bf = d - 1;
	g_bf_strand = gene_strand[g_bf];

	if(breakIT) {
		# previous IT on same strand
		it_bf = t-1; 
		while(it_strand[it_bf] && it_strand[it_bf] != g_bf_strand) { it_bf--; }
		if(!it_strand[it_bf]) { it_end[it_bf] = 0; }
	}

	# Initializing counter for nb of gene(s) before the IT on same strand
	nb_g_bf = 0;

	# Get info on the before series of genes
	while(1) {
		nb_g_bf++;
		bf_series[nb_g_bf] = write_gene_tag(g_bf);
		g_bf--;
		if(breakIT && gene_start[g_bf] <= it_end[it_bf]) {
			bf_break = "cause=IT;id=" it_id[it_bf] ";start=" it_start[it_bf] ";end=" it_end[it_bf];  
			break;
		}
		if(!gene_strand[g_bf]) {
			bf_break = "cause=First_gene_from_ori";
			break;
		}
		if(gene_strand[g_bf] != g_bf_strand ) {
			bf_break = "cause=Opposite_gene;" write_gene_tag(g_bf);
			break;
		}
	}

	# Concatenate the list of genes into one tag
	bf_series_tag = bf_series[nb_g_bf];
	for(i = nb_g_bf-1; i >= 1; i--) {
		bf_series_tag = bf_series_tag "|" bf_series[i];
	}
	delete bf_series;

	## AFTER SERIES ##

	g_af = d;
	g_af_strand = gene_strand[g_af];
	
	if(breakIT) {
		it_af = t+1; 
		while(it_strand[it_af] && it_strand[it_af] != g_af_strand) { it_af++; }
		if(!it_strand[it_af]) { it_start[it_af] = 10^9; }
	}

	nb_g_af = 0;
	while(1) {
		nb_g_af++;
		af_series[nb_g_af] = write_gene_tag(g_af);
		g_af++;
		if(breakIT && gene_end[g_af] >= it_start[it_af]) {
			af_break = "cause=IT;id=" it_id[it_af] ";start=" it_start[it_af] ";end=" it_end[it_af];
			break;
		}
		if(!gene_strand[g_af]) {
			bf_break = "cause=Last_gene_from_ori";
			break;
		}
		if(gene_strand[g_af] != g_af_strand) {
			af_break = "cause=Opposite_gene;" write_gene_tag(g_af);
			break;
		}

	}
	af_series_tag = af_series[1];
	for(i = 2; i <= nb_g_af; i++) {
		af_series_tag = af_series_tag "|" af_series[i];
	}
	delete af_series;

	genomic_landscape = g_bf_strand \
		"\t" nb_g_bf \
		"\t" bf_series_tag \
		"\t" bf_break \
		"\t" g_af_strand \
		"\t" nb_g_af \
		"\t" af_series_tag \
		"\t" af_break;

	return genomic_landscape;
}

BEGIN {
	FS = "\t";
	file_idx = 0;

	# Arg 'breakIT' is sent to this script
	# It is a boolean specifying whether the series of genes
	# can be broken by another IT or not

	# gene counter
	g = 0;
	# terminator counter
	t = 0;
	# counter for index of gene dictionary;
	#     meant to avoid redundant reading;
	d = 0;
}

FNR == 1 {
	file_idx++;
}

# The input file 1 is the annotation of the genome
# Lines starting with a comment are ignored
file_idx == 1 && FNR > 1 {
	gene_start[g] = $1;
	gene_end[g] = $2;
	gene_strand[g] = $3;
	gene_name[g] = $4;
	g++;

}

# The input file 2 is the list of terminators
file_idx == 2 && FNR == 1 {
	if(NF == 7) {
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
file_idx == 2 && FNR > 1 && NF == 7 {
	prefix[t] = $0;
	it_id[t] = $1;
	it_start[t] = $2;
	it_end[t] = $3;
	it_strand[t] = $4;
	t++;
}

# Otherwise the line is retrieved as such
file_idx == 2 && FNR > 1 && NF > 7 {
	already_line[t] = $0; t++;
}

END {
	printf("%s\n", header);
	for(i = 0; i < t; i++) {
		if(alread_line[i]) {
			printf("%s\n", already_line[i]);
		} else {
			# The fct get_genomic_landscape will increment 
			# the dictionary index 'd'
			suffix = get_genomic_landscape(i, d, breakIT);
			line = prefix[i] "\t" suffix;
			printf("%s\n", line);
		}
	}
}
