/*
THIS SCRIPT TAKES 2 INPUT FILES:
- FILE 1 MUST CONTAIN A LIST OF REFERENCE SEQUENCES
- FILE 2 MUST CONTAIN A LIST OF SEQUENCES TO BE COMPARED 
	 WITH EACH AND EVERY REFERENCE SEQUENCE

THE OUTPUT IS A TABLE OF HOMOLOGY
*/

//#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include "load_genes.h"
#include "edlib.h"


struct gene
{
  char *id;
  char *seq;
  struct gene *next;
};

const int nbAmino = 24;
struct blossum
{
  int scoreMatrix[nbAmino][nbAmino];
  int minScore;
  int maxScore;
};

struct fogsaaOut {
  int score;
  char *algA;
  char *algB;
};

static void
usage(char *argv0)
{
  fprintf(stderr, "%s expects three arguments:\n", argv0);
  fprintf(stderr,
	 "%s [ -r filename ] [ -c filename ] [ -o filename ]\n", argv0);
  fprintf(stderr,
	  "  -r: specify the file containing the reference list of sequences.\n");
  fprintf(stderr,
	  "  -c: specify the file containing the list of sequences to be compared with the reference list.\n");
  fprintf(stderr,
	  "  -o: specify the output file corresponding to the table of homology.\n");
}

void
compare(struct gene *refList, struct gene *vsList)
{
  char *refSeq, *vsSeq;
  int k = 0, refLen, vsLen, diffLen;
  struct gene *headVsList = vsList;
  
  while(refList != NULL) {
    
    refSeq = refList -> seq;
    refLen = strlen(refSeq);
    
    vsList = headVsList;
    while(vsList != NULL) {

      vsSeq = vsList -> seq;
      vsLen = strlen(vsSeq);
      
      diffLen = abs(refLen - vsLen);
      if(diffLen < (refLen / 4)) {
	EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen, edlibDefaultAlignConfig());
	printf("edit_distance(%s, %s) = %d\n", refList -> id, vsList -> id, result.editDistance);
	edlibFreeAlignResult(result);
      }

      vsList = vsList -> next;
    }
    refList = refList -> next;
  }
}

int
main(int argc, char **argv)
{
  int i;
  char *refFile, *vsFile, *outFile;
  
  struct gene *refListGenes, *vsListGenes;

  
  /* This scripts requires to be called with 7 args exactly:
  (1 is the script itself, 3 options and their 3 arguments) */
  if(argc != 7) {
    usage(argv[0]);
    exit(0);
  }
  i = 1;
  while(i < argc) {
    /* This script requires options ("-")" */
    if(argv[i][0] != '-') {
      usage(argv[0]);
      exit(1);
    }
    /* All options requires to specify an argument string (filename) */ 
    i++;
    if(i >= argc) {
      usage(argv[0]);
      exit(2);
    }
    /* Only 3 letters as options are allowed */
    switch(argv[i-1][1]) {
    	case 'r':
	  refFile = argv[i];
	  break;
    	case 'c':
	  vsFile = argv[i];
	  break;
    	case 'o':
	  outFile = argv[i];
	  break;
    	default :
	  usage(argv[0]);
	  exit(3); 
    }
    i++;
  }

  refListGenes = loadGenes(refFile);
  vsListGenes = loadGenes(vsFile);

  compare(refListGenes, vsListGenes);

  return 0;
}

				







