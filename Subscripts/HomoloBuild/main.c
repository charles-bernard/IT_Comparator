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
  fprintf(stderr,
	  "  -m: method of alignment (NW, HW, SHW).\n");
  fprintf(stderr,
	  "  -w: write pretty alignment of best match.\n");
  fprintf(stderr,
	  "  -a: avoid Redundant Comparison with already assigned and certained homologuous.\n");
}

/* This function is there in case one would like to print the alignment
   btw two homologuous proteins */
void
writePrettyAlignment(char *queryId, char* query, int queryLength,
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

/* if several vsSeqs happen to have the same best Score,
   this function is there to concatenate their characteristics */
char *
concat(const char *s1, const char *s2)
{
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);
    char *result = (char *)malloc(sizeof(char) * (len1+len2+2));//+2 for sep and zero-terminator
    memcpy(result, s1, len1);
    result[len1]=';';
    memcpy(result+len1+1, s2, len2+1);//+1 to copy the null-terminator
    return result;
}

/* This is the core function */
void
  compare(EdlibAlignMode mode, char *filename,
	  struct gene *refList, struct gene *vsList,
	  int avoidRedundance, int writeAlignment)
{
  int g = 1; int ng = lenList(refList);
  char *refSeq, *vsSeq, *bestSeq;
  char *refId, *vsId, *bestId = ".";
  int refLen, vsLen, bestLen, diffLen;
  
  struct gene *headVsList = vsList;
  struct gene *oldGene, *beforeBestGene, *bestGene;
  
  int score, bestScore;
  float scoreLenRatio;
  int bestPosition;
  int bestAlignmentLength;
  unsigned char *bestAlignment = (unsigned char *)malloc(sizeof(unsigned char));
  char *cigar = (char *)malloc(sizeof(char));

  FILE *f = fopen(filename, "w");
  if(f == NULL) {
    printf("Impossible to open \"%s\"\n", filename);
    exit(5);
  }
  fprintf(f, "refId\tvsId\tedlibAlignmentScore\tCigarAlignment\n");
  
  while(refList != NULL) {

    refId = refList -> id;
    refSeq = refList -> seq;
    refLen = strlen(refSeq);

    bestScore = refLen * 2;
    
    vsList = headVsList;
    if(avoidRedundance == 1) {
      oldGene =  beforeBestGene = vsList;
    }
          
    while(vsList != NULL) {

      vsId = vsList -> id;
      vsSeq = vsList -> seq;
      vsLen = strlen(vsSeq);
      
      diffLen = abs(refLen - vsLen);
      if(diffLen < (refLen / 5)) {

	/* Here is just a reminder of how to run edlib by default: */
	//EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen, edlibDefaultAlignConfig());

	// NB the major difference with main.c (for genes) is the "HW" infix EDLIB_MODE
	EdlibAlignResult result = edlibAlign(refSeq, refLen, vsSeq, vsLen,
					     edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH));
	score = result.editDistance;
	
	if(result.alignment && score <= bestScore) {
	  
	  /* Store variables relative to this new best target */
	  if(score < bestScore) {
	    bestId = vsId;
	    cigar = edlibAlignmentToCigar(result.alignment,
					  result.alignmentLength,
					  EDLIB_CIGAR_EXTENDED);
	  } else if(score == bestScore) {
	    bestId = concat(bestId, vsId);
	    cigar = concat(cigar,
			   edlibAlignmentToCigar(result.alignment,
						 result.alignmentLength,
						 EDLIB_CIGAR_EXTENDED));
	  }

	  if(avoidRedundance == 1) {
	    bestGene = vsList;
	    beforeBestGene = oldGene;
	  }
	  
	  if(writeAlignment == 1) {
	    free(bestAlignment);
	    bestSeq = vsSeq;
	    bestLen = vsLen;
	    bestPosition = *result.endLocations;
	    bestAlignmentLength = result.alignmentLength;
	    bestAlignment = (unsigned char *)malloc(sizeof(unsigned char) * (bestAlignmentLength));
	    memcpy(bestAlignment, result.alignment, (bestAlignmentLength));
	  }
	  
	  bestScore = score;
	}
	edlibFreeAlignResult(result);
      }

      if(avoidRedundance == 1) {
	oldGene = vsList;
      }
      vsList = vsList -> next;
    }
    
    /* The task might be quite long, so a counter is necessary */ 
    printf("\rTreated Query Seq n°%d/%d", g, ng);
    fflush(stdout);
    
    scoreLenRatio = float(bestScore) / float(refLen);
    if(scoreLenRatio < 1) { //0.5 is found to be a good threshold
      fprintf(f, "%s\t%s\t%.2f\t%s\n", refId, bestId, scoreLenRatio, cigar);

      if(writeAlignment == 1) {
	writePrettyAlignment(refId, refSeq, refLen, bestId, bestSeq, bestLen,
		       bestPosition, bestAlignment, bestAlignmentLength);
      }
      
      if(avoidRedundance == 1) {
	if(scoreLenRatio < 0.1 && headVsList != bestGene) {
	  beforeBestGene -> next = bestGene -> next;
	  free(bestGene);
	}
      }

    } else
      fprintf(f, "%s\t.\t.\t.\n", refId);
    
    refList = refList -> next;
    g++;
  }
  free(cigar);
  free(bestAlignment);
  printf("\n");
  fclose(f);
}

int
main(int argc, char **argv)
{
  int i;
  char *refFile, *vsFile, *outFile;
  struct gene *refListGenes, *vsListGenes;

  char *method = "NW";
  EdlibAlignMode mode;

  int writeAlignment = 0, avoidRedundance = 0;
  
  /* This scripts requires to be called with 7 args exactly:
  (1 is the script itself, 3 options and their 3 arguments) */
  if(argc < 7) {
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
    i++;
    /* All options requires to specify an argument string (filename)  
    i++;
    if(i >= argc) {
      usage(argv[0]);
      exit(2);
    }
    */
    /* Only 6letters as options are allowed */
    switch(argv[i-1][1]) {
    case 'r':
      refFile = argv[i];
      i++;
      break;
    case 'c':
      vsFile = argv[i];
      i++;
      break;
    case 'o':
      outFile = argv[i];
      i++;
      break;
    case 'm':
      method = argv[i];
      i++;
      break;
    case 'w':
      writeAlignment = 1;
      break;
    case 'a':
      avoidRedundance = 1;
      break;
    default :
      usage(argv[0]);
      exit(3); 
    }
  }

  if(strcmp(method, "NW") == 0) {
    mode = EDLIB_MODE_NW;
  } else if(strcmp(method, "HW") == 0) {
    mode = EDLIB_MODE_HW;
  } else if(strcmp(method, "SHW") == 0) {
    mode = EDLIB_MODE_SHW;
  } else {
    fprintf(stderr, "Wrong method of alignment!\n");
    usage(argv[0]);
    exit(4);
  }
      
  refListGenes = loadGenes(refFile);
  vsListGenes = loadGenes(vsFile);

  compare(mode, outFile, refListGenes, vsListGenes, avoidRedundance, writeAlignment);

  freeList(refListGenes);
  freeList(vsListGenes);

  return 0;
}
