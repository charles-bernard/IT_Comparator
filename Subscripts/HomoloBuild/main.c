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
printAlignment(const char* query, const int queryLength,
	       const char* target, const int targetLength,
	       EdlibAlignResult result)
{
  int position = *result.endLocations;
  unsigned char *alignment = result.alignment;
  int alignmentLength = result.alignmentLength;
  
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
compare(struct gene *refList, struct gene *vsList)
{
  char *refSeq, *vsSeq, *bestSeq;
  int k = 0, refLen, vsLen, bestLen, diffLen;
  struct gene *headVsList = vsList;

  int score, bestScore;
  char *bestId;
  //EdlibAlignResult bestResult;
  
  while(refList != NULL) {
    
    refSeq = refList -> seq;
    refLen = strlen(refSeq);
    
    vsList = headVsList;
    while(vsList != NULL) {

      bestScore = -1000;

      vsSeq = vsList -> seq;
      vsLen = strlen(vsSeq);
      
      diffLen = abs(refLen - vsLen);
      if(diffLen < (refLen / 4)) {
	//EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen, edlibDefaultAlignConfig());
	EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH));
	score = result.editDistance;
	printf("%d\n", score);
	printAlignment(refSeq, refLen, vsSeq, vsLen, result);
	edlibFreeAlignResult(result);

	/*
	if(score > bestScore) {
	  bestId = vsList -> id;
	  bestScore = score;
	  //edlibFreeAlignResult(bestResult);
	  bestSeq = vsSeq;
	  bestLen = vsLen;
	  //EdlibAlignResult bestResult = result;
	  edlibFreeAlignResult(result);
	} else {
	  edlibFreeAlignResult(result);
	}
	*/
	
      }
      vsList = vsList -> next;
    }

    /*
    printf("Best Score:%d\n", bestScore);
    printf("RefId: %s;BestId: %s\n", refList -> id, bestId);
    printf("RefSeq:%s\n", refSeq);
    printf("BestSeq:%s\n", bestSeq);
    */
    
    //printAlignment(refSeq, refLen, bestSeq, bestLen, bestResult);
    //edlibFreeAlignResult(bestResult);
    
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

				







