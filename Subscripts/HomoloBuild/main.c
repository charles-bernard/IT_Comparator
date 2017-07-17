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
#include "list_genes.h"
#include "edlib.h"

struct gene
{
  char *id;
  char *seq;
  struct gene *next;
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

/* This function is there in case one would like to print the alignment
   btw two homologuous proteins */
void
writeAlignment(char *queryId, char* query, int queryLength,
	       char *targetId, char* target, int targetLength,
	       int position, unsigned char *alignment,
	       int alignmentLength)
{
  printf("(Q: %s) VS (T: %s)\n\n", queryId, targetId);

  int tIdx = -1;
  int qIdx = -1;
  for (int start = 0; start < alignmentLength; start += 50) {
    // target
    printf("T: ");
    int startTIdx;
    for (int j = start; j < start + 50 && j < alignmentLength; j++) {
      if (alignment[j] == EDLIB_EDOP_INSERT)
	printf("-");
      else
	printf("%c", target[++tIdx]);
      if (j == start)
	startTIdx = tIdx;
    }
    if(startTIdx > 0)
      printf(" (%d - %d)\n", startTIdx, tIdx);
    else
      printf(" (%d - %d)\n", 0, tIdx);

    // match / mismatch
    printf("   ");
    for (int j = start; j < start + 50 && j < alignmentLength; j++) {
      printf(alignment[j] == EDLIB_EDOP_MATCH ? "|" : " ");
    }
    printf("\n");

    // query
    printf("Q: ");
    int startQIdx = qIdx;
    for (int j = start; j < start + 50 && j < alignmentLength; j++) {
      if (alignment[j] == EDLIB_EDOP_DELETE)
	printf("-");
      else
	printf("%c", query[++qIdx]);
      if (j == start)
	startQIdx = qIdx;
    }
    if(startQIdx > 0)
      printf(" (%d - %d)\n", startQIdx, qIdx);
    else
      printf(" (%d - %d)\n", 0, qIdx);
  }
}

void
compare(char *filename, struct gene *refList, struct gene *vsList, int wa)
{
  int g = 1; int ng = lenList(refList);
  char *refSeq, *vsSeq, *bestSeq;
  char *refId, *vsId, *bestId = ".";
  int refLen, vsLen, bestLen, diffLen;
  
  struct gene *headVsList = vsList;
  struct gene *oldGene = vsList, struct gene *previousBestGene = vsList;
  struct gene *bestGene;
  
  int score, bestScore;
  float scoreLenRatio;
  int bestPosition;
  int bestAlignmentLength;
  unsigned char *bestAlignment = (unsigned char *)malloc(sizeof(unsigned char));

  FILE *f = fopen(filename, "w");
  fprintf(f, "refId\tvsId\tedlibAlignmentScore\n");
  
  while(refList != NULL) {

    refId = refList -> id;
    refSeq = refList -> seq;
    refLen = strlen(refSeq);

    bestScore = refLen * 2;
    
    vsList = headVsList;
          
    while(vsList != NULL) {

      vsId = vsList -> id;
      vsSeq = vsList -> seq;
      vsLen = strlen(vsSeq);
      
      diffLen = abs(refLen - vsLen);
      if(diffLen < (refLen / 5)) {

	/* Here is just a reminder of how to run edlib by default: */
	//EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen, edlibDefaultAlignConfig());
	
	EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen,
					     edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
	score = result.editDistance;	  
	if(result.alignment && score < bestScore) {
	  bestScore = score;
	  bestId = vsId;
	  bestGene = vsList;
	  previousBestGene = oldGene;
	  if(wa == 1) {
	    free(bestAlignment);
	    bestSeq = vsSeq;
	    bestLen = vsLen;
	    bestPosition = *result.endLocations;
	    bestAlignmentLength = result.alignmentLength;
	    bestAlignment = (unsigned char *)malloc(sizeof(unsigned char) * (bestAlignmentLength));
	    memcpy(bestAlignment, result.alignment, (bestAlignmentLength));
	  }
	}
	edlibFreeAlignResult(result);
      }
      oldGene = vsList;
      vsList = vsList -> next;
    }
    /* The task might be quite long, so a counter is necessary */ 
    printf("\rtreated gene n°%d/%d", g, ng);
    fflush(stdout);
    
    scoreLenRatio = float(bestScore) / float(refLen);
    if(scoreLenRatio < 0.5) {
      if(wa == 1) {
	writeAlignment(refId, refSeq, refLen, bestId, bestSeq, bestLen,
		       bestPosition, bestAlignment, bestAlignmentLength);
      }
      fprintf(f, "%s\t%s\t%.2f\n", refId, bestId, scoreLenRatio);
      
      //prevent redundant comparison with already assigned and certain homologuous
      if(scoreLenRatio < 0.1 && headVsList != bestGene) {
	previousBestGene -> next = bestGene -> next;
	free(bestGene);
      }

    } else
      fprintf(f, "%s\t.\t.\n", refId);
    refList = refList -> next;
    g++;
  }
  free(bestAlignment);
  fclose(f);
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

  compare(outFile, refListGenes, vsListGenes, 0);

  freeList(refListGenes);
  freeList(vsListGenes);

  return 0;
}

				







