#ifndef LOAD_GENES_H_
#define LOAD_GENES_H_

void testFile(char *filename);
struct fields *getFields(char *line);
struct gene *appendList(struct gene *list, char *id, char *seq);
void printList(struct gene *list);
void freeList(struct gene *list);
struct gene *loadGenes(char *filename);

#endif // LOAD_GENES_H_
