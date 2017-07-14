/*
THIS SCRIPT TAKES 2 INPUT FILES:
- FILE 1 MUST CONTAIN A LIST OF REFERENCE SEQUENCES
- FILE 2 MUST CONTAIN A LIST OF SEQUENCES TO BE COMPARED 
	 WITH EACH AND EVERY REFERENCE SEQUENCE

THE OUTPUT IS A TABLE OF HOMOLOGY
*/

//#include <iostream.h>
#include <stdio.h>;
#include <string.h>;
#include <stdlib.h>;
#include <sys/time.h>;
#include <unistd.h>;


struct fields
{
  char *id;
  char *seq;
};


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


static void
testFile(char *filename)
{
  FILE *f;
  f = fopen(filename, "r");
  if(f == 0) {
    fprintf(stderr,
	    "Error: Impossible to open \"%s\"\n", filename);
    exit(4);
  }
  fclose(f);
}


struct fields *
getFields(char *line)
{
  char *sep = "\t\n";
  char *p;
  int k;
  struct fields *fields = (struct fields *)malloc(sizeof(struct fields));
  
  p = strtok(line, sep);
  k = 1;
  while(p != NULL) {
    if(k == 1) {
      fields -> id = strdup(p);
      k = 2;
    } else if(k == 2) {
      fields -> seq = strdup(p);
      k = 1;
    }
    p = strtok(NULL, sep);
  }
  return fields;
}
  

struct gene *
append(struct gene *list, char *id, char *seq)
{
  struct gene *newGene = (struct gene *)malloc(sizeof(struct gene));
  list -> next = newGene;
  newGene -> id = id;
  newGene -> seq = seq;
  newGene -> next = NULL;
  return newGene;
}


void
printList(struct gene *list)
{
  if(list == NULL)
    printf("Empty list of genes\n");
  else {
    printf("%s:%s\n", list -> id, list -> seq);
    if(list -> next != NULL)
      printList(list -> next);
  }
}


void
freeList(struct gene *list)
{
  if(list -> next != NULL) {
    freeList(list -> next);
    free(list);
  }
}


struct gene *
parseFile(char *filename)
{
  FILE *f = fopen(filename, "r");
  
  const size_t bufSize = 4096; // This should be enough for protein sequences
  char buf[bufSize];
  struct fields *fields;
  
  int isHead = 1;
  struct gene *head = (struct gene *)malloc(sizeof(struct gene));
  struct gene *list = (struct gene *)malloc(sizeof(struct gene));

  while(fgets(buf, bufSize, f)) {
    fields = getFields(buf);
    list = append(list, fields -> id, fields -> seq);
    if(isHead == 1) {
      head = list;
      isHead = 0;
    }
  }
  fclose(f);
  return(head);
}


void
compare(struct gene *refList, struct gene *vsList)
{
  struct gene *headVsList = vsList;
  while(refList != NULL) {
    vsList = headVsList;
    while(vsList != NULL) {
      printf("%s:%s///%s:%s\n",
	     refList -> id, refList -> seq,
	     vsList -> id, vsList -> seq);
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

  testFile(refFile);
  testFile(vsFile);

  parseFile(refFile);
  //parseFile2(vsFile);
  refListGenes = parseFile(refFile);
  //vsListGenes = parseFile(vsFile);

  //compare(refListGenes, vsListGenes);
  
  printList(refListGenes);
  //printList(vsListGenes);
  
  return 0;
}

				







