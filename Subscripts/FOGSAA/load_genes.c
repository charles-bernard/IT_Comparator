#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "load_genes.h"

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

void
testFile(char *filename)
{
  FILE *f;
  f = fopen(filename, "r");
  if(f == 0) {
    fprintf(stderr, "Error: Impossible to open \"%s\"\n", filename);
    exit(4);
  }
  fclose(f);
}

struct fields *
getFields(char *line)
{
  char sep[4] = "\t\n";
  char *p;
  int k = 1;
  struct fields *fields = (struct fields *)malloc(sizeof(struct fields));
  
  p = strtok(line, sep);
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
appendList(struct gene *list, char *id, char *seq)
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
  }
  free(list);
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
    list = appendList(list, fields -> id, fields -> seq);
    if(isHead == 1) {
      head = list;
      isHead = 0;
    }
  }
  fclose(f);
  return(head);
}
