//#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include "list_genes.h"
#include "edlib.h"


int
main(int argc, char **argv)
{
  char *query;
  char *target;

  int queryLen;
  int targetLen;
  
  EdlibAlignMode mode = EDLIB_MODE_HW;

  float scoreLenRatio = 1.0;
  
  query = argv[1];
  target = argv[2];

  queryLen = strlen(query);
  targetLen = strlen(target);
  
  EdlibAlignResult result = edlibAlign(query, queryLen, target, targetLen,
				       edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH));

  scoreLenRatio = float(result.editDistance) / float(queryLen);
  printf("%.2f", scoreLenRatio);
  
  edlibFreeAlignResult(result);

  return 0;
}
