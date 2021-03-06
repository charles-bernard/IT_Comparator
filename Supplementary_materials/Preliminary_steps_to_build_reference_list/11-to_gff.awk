#/usr/bin/awk

BEGIN {
	FS = "\t";
	source = "IT_Miner";
	# the variables 'organism' is sent to this script;
	if(!organism) {
		organism = "whole genome";
	};
	feature = "terminator";
	strength = ".";

	header = "##<organism>" \
		"\t<source>" \
		"\t<feature>" \
		"\t<start>" \
		"\t<end>" \
		"\t<strand>" \
		"\t<RNIE_bit_score>" \
		"\t<.>" \
		"\t<comment>"
	#print header; 
}

NR > 1 {
	start = $1;
	end = $2;
	strand = $3;
	score = ".";
	sequence = $4;
	energy = $6;
	structure = $5;

	if(strand == "+") {
		id = "IT-" end;
	} else {
		id = "IT-" start;
	}

	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t", 
		organism,
		source,
		feature,
		start,
		end,
		score,
		strand)
	printf("ID=%s;sequence=%s;pred_structure=%s;pred_free_energy=%s;pred_strength=%s\n",
		id,
		sequence,
		structure,
		energy,
		strength)
}
