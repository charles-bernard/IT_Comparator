#!/usr/bin/awk

BEGIN {
	FS = "\t";
	strand["forward"] = "+";
	strand["reverse"] = "-";
	species = "ECK12";
	source = "RegulonDB";
	feature = "IT";
}

!/^#/ && $7 == "rho-independent" {
	if($4=="both") {
		printf("%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t.\n", \
			species, source, feature, \
			$2, $3, "+");
		printf("%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t.\n", \
			species, source, feature, \
			$2, $3, "-");
	} else {
		printf("%s\t%s\t%s\t%s\t%s\t.\t%s\t.\t.\n", \
			species, source, feature,
			$2, $3, strand[$4]);
	}
}