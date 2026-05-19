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
  char qfilename[20];
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
