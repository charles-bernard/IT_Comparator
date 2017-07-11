#!/usr/bin/awk

function get_tags(attribute) {
	match(attribute, /ID=[^;]*;/, id_tag);
	match(attribute, /sequence=[^;]*;/, sequence_tag);
	match(attribute, /pred_structure=[^;]*;/, structure_tag);
	match(attribute, /pred_free_energy=[^;]*;/, energy_tag);

	gsub(/^.*=/, "", id_tag[0]); gsub(/;$/, "", id_tag[0]);
	gsub(/^.*=/, "", sequence_tag[0]); gsub(/;$/, "", sequence_tag[0]);
	gsub(/^.*=/, "", structure_tag[0]); gsub(/;$/, "", structure_tag[0]);
	gsub(/^.*=/, "", energy_tag[0]); gsub(/;$/, "", energy_tag[0]);

	return id_tag[0] ";" sequence_tag[0] ";" structure_tag[0] ";" energy_tag[0];
}

BEGIN {
	FS = "\t";
	header = "ID" \
		"\tStart" \
		"\tEnd" \
		"\tInitially Detect Strand" \
		"\tPrimary Sequence" \
		"\tSecondary Structure (Dot-Bracket Notation)" \
		"\tFree Energy (kcal/mol)"
	printf("%s\n", header);
}

!/^#.*$/ {
	start = $4;
	end = $5;
	strand = $7;

	attribute = $9;
	attribute_tags = get_tags(attribute);
	split(attribute_tags, fields, ";");

	id = fields[1];
	sequence = fields[2];
	structure = fields[3];
	energy = fields[4];

	printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
		id, start, end, strand, sequence, structure, energy);
}

