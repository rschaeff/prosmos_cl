#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include "h_bond.h"
#include "Fresidue.h"
#include "control.h"
#include "element.h"
#include "mpi.h"
using namespace std;
int main(int argc, char* argv[])
{
  int i;
  int debugk = 0;
  char tokenstring[50];
  FILE *flist,*a,*hander;
  char fileName[50];
  vector<sheet> sheetset;
  bool alternative = false;
  int numprocs,myrank;
  control con;
  vector<element> eleInOneStru;
  char fileName2[300];
  char dirname[500];
  //char tokenstring[200];

  MPI_Init(&argc ,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  if(argc<4 ||argc>5 || strcmp(argv[1],"-o")!=0 && strcmp(argv[1],"-f")!=0 && strcmp(argv[1],"-d")!=0&&
        strcmp(argv[1],"-os")!=0&&strcmp(argv[1],"-fs")!=0&&strcmp(argv[1],"-ds")!=0)
  {
     if(myrank == 0){
     cout<<"we have three format of this command line"<<endl;
     cout<<"the first format"<<endl;
     cout<<"this format of the command line will generate all the interaction matrix for all pdb files"<<endl;
     cout<<"the first paremeter is the name of the object file"<<endl;
     cout<<"the second parameter [-fs] is the file option"<<endl;
     cout<<"the program will generate the sheet information by reading PALSSE files"<<endl;
     cout<<"the third parameter is the name of the total pdb file name list file and path name"<<endl;
     cout<<"the fourth parameter is only the path name of PALSSE outfile"<<endl;
     cout<<"the fifth parameter is the path name and file name of your output file"<<endl;
     cout<<"for example :a.out -fs matrixsample1 ../ss-vector/ ./matrixOutput"<<endl;
     cout<<"the second format : "<<endl;
     cout<<"this command will generate all the interaction matrix for all pdb files in one directory"<<endl;
     cout<<"the first paremeter is the name of the object file"<<endl;
     cout<<"the second parameter [-ds] is the direstory option"<<endl;
     cout<<"the program will generate the sheet information by reading PALSSE files"<<endl;
     cout<<"the third parameter is only the path name of all pdb files"<<endl;
     cout<<"the fouth parameter is the path name and file name of your output file"<<endl;
     cout<<"for example :a.out -ds  ../ss-vector/ ./matrixOutput"<<endl; 
     cout<<"this is the third format :"<<endl;
     cout<<"this is format of the five parameter"<<endl;
     cout<<"this format of the command line will generate one interaction matrix for one pdb file"<<endl;
     cout<<"the first paremeter is the name of the object file"<<endl;
     cout<<"the second parameter is the name of pdb id"<<endl;
     cout<<"the third parameter is an fixed option -os"<<endl;
     cout<<"the program will generate the sheet information by reading PALSSE file"<<endl;
     cout<<"the fouth parameter is only the path name of PALSSE outfile"<<endl;
     cout<<"the fifth parameter is the path name and file name of your output file"<<endl;
     cout<<"for example :a.out -os 1b5s.ssd ../ss-vector/ ./singlepdb"<<endl;
     }
     alternative = true;
     MPI_Finalize();
     return 0;
  }
  if(argc == 5 || argc == 4){
  if(strcmp(argv[1],"-d")==0||strcmp(argv[1],"-fs")==0||strcmp(argv[1],"-ds")==0)
  {
     con.setoption(argv[1]);
     if(myrank == 0)
       con.manager(numprocs,argv[1],argv[2]);
     else
     {
        if(argc == 5)
          con.worker(argv[3],argv[4]);
        else
          con.worker(argv[2],argv[3]);
     } 
  }
  if(strcmp(argv[1],"-o")==0||strcmp(argv[1],"-os")==0) 
  {
      if(myrank == 0)
      {
          eleInOneStru.clear();
          element myelement;
          a = fopen(argv[4] ,"w");
          if(a == NULL)
          {
            cout<<argv[4]<<" can't open "<<endl;
            MPI_Finalize();
            return 0;
          }
          con.setoption(argv[1]);
          strcat(argv[3],argv[2]);
          strcpy(fileName,argv[3]);
          strcpy(tokenstring ,argv[2]); 
          cout<<"the file name is "<<fileName<<endl;
          con.prepareIndex(fileName);
          con.readIndrefile1(fileName , eleInOneStru);
          cout<<"this out of read file"<<endl;
          //cout<<"the eleInOneStru size is "<<eleInOneStru.size()<<endl;
          con.Elmentdelete(eleInOneStru);
          con.prodvector(eleInOneStru);
          con.proInteractionMatr(eleInOneStru,a,tokenstring);
          fclose(a);
          MPI_Finalize();
          return 0;
      }
   }
   if(strcmp(argv[1],"-f") == 0)
   {
      if(myrank == 0)
      {
        con.setoption(argv[1]);
        strcpy(fileName2,argv[2]);
        hander = fopen(argv[4],"w");
        if(a == NULL)
        {
            cout<<argv[4]<<" can't open "<<endl;
            MPI_Finalize();
            return 0;
        }
        ifstream inputfile(fileName2,ios::in);
        if(!inputfile)
        {
           cout<<"the file "<<fileName2<<" can't open"<<endl;
           MPI_Finalize();
           exit(0);
        }
        else
        {
           while (inputfile.getline(tokenstring ,20,'\n'))
           {
              eleInOneStru.clear();  
              strcpy(dirname,argv[3]);
              strcat(dirname,tokenstring);
              cout<<"the dirname is "<<dirname<<endl;
              con.prepareIndex(fileName);
              con.readIndrefile1(dirname , eleInOneStru);
              con.Elmentdelete(eleInOneStru);
              con.prodvector(eleInOneStru);
              con.proInteractionMatr(eleInOneStru,hander,tokenstring);
           }
        }
     }
   }    
 }
 MPI_Finalize();
 return 0;
}
