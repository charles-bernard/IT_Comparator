#include <stdio.h>
#include <string.h> 
#include <stdlib.h>
#include "load_amino_and_blossum.h"

const int nbAmino = 24;

struct blossum
{
  int scoreMatrix[nbAmino][nbAmino];
  int minScore;
  int maxScore;
};
  
char *
loadAmino(char *filename)
{
  FILE *f = fopen(filename, "r");
  char *amino = (char *)malloc(sizeof(char) * (nbAmino + 1));
  char c;
  int i;
  for(i = 0; i < nbAmino; i++) {
    fscanf(f, "%c ", &c);
    amino[c-'A'] = i;
  }
  fclose(f);
  return(amino);
}

struct blossum *
loadBlossum(char *filename)
{
  FILE *f;
  f = fopen("fscore.txt","r");
  int scoreMatrix[nbAmino][nbAmino];
  int minScore, maxScore;
  struct blossum *blossum = (struct blossum *)malloc(sizeof(struct blossum));
  
  fscanf(f, "%d", &scoreMatrix[0][0]);
  minScore = scoreMatrix[0][0];
  maxScore = scoreMatrix[0][0];
  
  fseek(f, 0L, SEEK_SET);
  for(int i = 0; i < nbAmino; i++) {
    for(int j = 0; j <= i; j++) {
      fscanf(f, "%d", &scoreMatrix[i][j]);
      scoreMatrix[j][i] = scoreMatrix[i][j];
      if(scoreMatrix[i][j] > maxScore)
	maxScore = scoreMatrix[i][j];
      if(scoreMatrix[i][j] < minScore)
	minScore = scoreMatrix[i][j];
    }
  }
  blossum -> minScore = minScore;
  blossum -> maxScore = maxScore;
  memcpy(blossum -> scoreMatrix, scoreMatrix, nbAmino * nbAmino * sizeof(int));
  fclose(f);
  return(blossum);
}
