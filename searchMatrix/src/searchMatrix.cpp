#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <iostream>
#include <vector>
#include <queue>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#ifdef SILENT
#include <cstdio>
#endif
#include "elecol.h"
#include "searchControl.h"
#include "Fpass.h"
#include "external.h"
#include "MtrixElment.h"
#include "handness.h"
#include "sheet.h"
#include "require.h"
using namespace std;
int main(int argc , char *argv[])
{
#ifdef SILENT
  // searchControl.h has 200+ leftover debug cout/printf calls inside the
  // per-DB-entry match loop; on a 710k-entry DB that's tens of millions of
  // small writes per query. Redirect stdout once at startup -- both cout and
  // printf hit fd 1, so this catches them all without touching the call sites.
  // cerr (fd 2) is preserved so real errors still surface. usage/arg-error
  // messages above were also cout, but those exit(0) anyway and the caller
  // only inspects hit count + rc, so losing them is fine.
  FILE *silent_fp = freopen("/dev/null", "w", stdout);
  (void)silent_fp;
#endif
  bool judge1 = false;
  bool judge2 = false; 
  char **interActionM;
  char **quryMatrix;
  vector<elecol> intElecol;
  vector<elecol> qElecol;
  vector<handness> totalhandqury;
  vector<require> totalrequire;
  vector<sheet> totalsheet;
  elecol *eptr;
  char qfilename[4096];  /* was [20]: overflowed on query paths >19 chars (stack smash) */
  char quline[1000];
  char intMline[1000000];
  FILE *qufevt, *mtfevt;
  searchControl control;
  int i , j, qrow,mrow;
  int readLincon;
  vector<matrixElment> IntMnumEl;
  sheet *sheetptr;
  vector<fpass> totalpass;
  char pid[20];
  int argvint=0;
  char command[200];
  int numprocs;
  int myrank;
  if(argc<4)
  {
      cout<<"command line should have three parameters"<<endl;
      cout<<"the first parameter is name of the object file"<<endl;
      cout<<"the second parameter is the path and file name of qury file"<<endl;
      cout<<"the third parameter is the path and out put generate matrix file"<<endl;
      cout<<"the fouth parameter is the path name for the hit file "<<endl;
      exit(0);
  } 
  cout<<"this is step one : "<<endl;
  control.oneprocess(argv[1], argv[2] ,argv[3]);
  return 0;
}
