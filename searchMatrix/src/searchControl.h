#ifndef SEARCHCONTROL_H
#define SEARCHCONTROL_H
#include <queue>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#if defined (BSD) && !_POSIX_SOURCE
    #include <sys/dir.h>
#else
    #include <dirent.h>
#endif
#include <dirent.h>
#include "Fpass.h"
#include "MtrixElment.h"
#include "handness.h"
#include "sheet.h"
#include "require.h"
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;
class searchControl
{
  private :
          bool chainjudge;
          int countpid;
  public  :
       searchControl();
       void searchM(vector<elecol> a , char ** intM, int row1,vector<elecol>b,char** quM,int row2,vector<fpass>& totalpass
                    ,vector<matrixElment> &totalele);
       bool compareCol(vector<int> allRow,char **intM,int row1,char **quM,int row2,int nocol);
       void processMatrfile(char *filename,char **matrix);
       //this function is used to get the number of the element, so I can get the
       //the qury matrix of the element type and the corresponding colum
       int  quNumOfele(char *line, vector<elecol> &quelecol);
       //this function is used to get the qury matrix,which get the data from the 
       // the qury matrix file.
       void formqMatrix(char ** matrix,int row,ifstream &a,vector<handness> &hand ,vector<require> &tR);
       // this function I can get the number of the element in the intaction matrix
       // and the cordinates of the vector.
       // parameter a get from the main function:, oneline parameter is got from the
       // intactioMatrix
       // parameter pid is past from the main funtion is store for the structure id
       int intMnumofele(vector<matrixElment> &a, char *oneline,vector<elecol> &intEle,char *pid);
       //this function I will get the interaction matrix from yioutputmatrix, and fill the
       // the two dimension matrix interActionM.
       void getInterActionM(char ** intM,int mrow , char *oneline,vector<sheet> &totalsheet);
       //this function is used to judge if the two character is same or not,
       //because the qury matrix we have different definations about the letter.
       bool notequal(char a , char b);
       //this function is used to get the handness,give three elements to 
       //judge if the three element are right or left hand
       char chirality(matrixElment seg1,matrixElment seg2,matrixElment seg3);
       void selectMatrix(vector<matrixElment> &a,vector<handness> &hand,vector<require> &tR,vector<fpass>&allps,vector<sheet> &totalsheet,vector<elecol> &intEle,char **intM);
       bool sameOrdiffsheet(vector<sheet> &totalsheet,vector<int>eleid);
       void printOuptfile(vector<fpass> &totallpass,vector<matrixElment> &IntMnumEl,char *pid,char *path);
       void checkNumberLine(char *quline);
       bool formajob(char *filelist,char *job,int &size);
       void oneprocess(char *a1 , char *a2 ,char *a3);
       void setjudge(bool a);
       void checksheetH(char ** intactionM, int row,vector<sheet>&totalsheet,char *pid,vector<elecol> &intEle);
       int paraOrantisearch(vector<elecol> &intEle , int source ,int destinationtor, vector<int>&colum,vector<int> ssheet,
                            char **intM);
       void printgraph(vector<elecol> intEle ,int sou , int des ,vector<int> &pah);
};

searchControl::searchControl()
{
  chainjudge = true  ;
}
void searchControl::setjudge(bool a)
{
  chainjudge = a;
}

//this function is to check if the information in the sheet if correct or not, and also form the graph informatio
//to be used in the parallel or antiparallel search
void searchControl::checksheetH(char **intactionM,int row,vector<sheet>&totalsheet,char *pid,vector<elecol>&intEle)
{
   int numbersheet;
   int i,j,z,t,elecolindex;
   int rownum , colnum;
   numbersheet = totalsheet.size();
   char command[200];
   char sheetfilename[100];
   vector<int> sheetid;
   FILE *fvetc;
   bool fileopen = false;
   bool find = false;
   int w;
   bool judge = false;
   strcpy(sheetfilename ,"../sheetbug/");
   strcat(sheetfilename , pid);
   strcat(sheetfilename ,"sheet.bug");
   // checksheetH runs once PER DB RECORD, so a system("mkdir") here was a
   // fork+exec of mkdir for every one of ~5M records -- ~92% of per-record
   // runtime (12x slowdown). The directory only needs to exist once; make it
   // once per process.
   static bool sheetdir_made = false;
   if(!sheetdir_made)
   {
      strcpy(command, "mkdir  -p  ../sheetbug/" );
      system(command);
      sheetdir_made = true;
   }
   char sheetbugfilename[50];
   strcpy(sheetbugfilename , "../sheetbug/total.txt");
   fvetc = fopen(sheetbugfilename, "a");
   fileopen = true;
   if(fvetc == NULL)
   {
       cout<<"the file "<<sheetfilename<<" can't open "<<endl;
       exit(0);
   }
   for(i=0;i<numbersheet;i++)
   {
       sheetid.clear();
       sheetid = totalsheet[i].getSheet();
       if(sheetid.size()>=2)
       {
           for(j=0 ; j<sheetid.size() ;j++)
           {
               //this step I will find the proper element colum then add neighbor to it
               find = false;
               for(t=0;t<intEle.size();t++) 
               {
                   if(intEle[t].getcol() == sheetid[j])
                   {
                       elecolindex = t;
                       find = true;
                       break;
                   }
               }
               if(find == false)
               {
                   cout<<"in the checksheetH funtion, the intEle data or sheetid data is wrong"<<endl;
                   exit(0);
               }
               judge = false;
               for(z=0;z<sheetid.size();z++)
               {
                  if(intactionM[sheetid[j]][sheetid[z]] == 'c' || intactionM[sheetid[j]][sheetid[z]] == 't')
                  {
                     // intEle[sheetid[z]].setparalleloranti(intactionM[sheetid[j]][sheetid[z]]);
                      intEle[elecolindex].addneighbor(sheetid[z]);
                      judge = true;
                      //break;
                  }
               }
               /*
               if(judge == true)
               {
                  intEle[elecolindex].print();
                  exit(0);
               } 
               */
               if(judge == false)
               {
                  fprintf(fvetc," the no %s %d strand is not H with any other strand\n\n",pid,sheetid[j]+1);
               }
           }
           /*
           for(t=0;t<sheetid.size();t++)
           {
               if(intEle[sheetid[t]].getneighbor().size()==0)
               {
                   cout<<"bug in the sheet information"<<endl;
                   exit(0);
               }
           } 
          */ 
       } 
   }
    if(fileopen == true)
      fclose(fvetc);
} 

void searchControl::oneprocess(char *a1 , char *a2 ,char *a3)
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
  char qfilename[4096];  /* was [20]: overflowed on query paths >19 chars (stack smash) */
  char quline[1000];
  FILE *qufevt, *mtfevt;
  char intMline[1000000];
  //char intMline[100];
  int i , j, qrow,mrow;
  int readLincon;
  char logstring[5000];
  vector<matrixElment> IntMnumEl;
  sheet *sheetptr;
  vector<fpass> totalpass;
  char logfilename[100];
  FILE *logptr;
  char pid[20];
  int argvint=0;
  char command[200];
  char token[500];
  dirent *listptr;
  DIR   *dp;
  char dirname[500];
  string temp;
  string check;
  string tempstring;
  pid_t mypid;
  char command1[500];
  char command2[500];
  int indexc = 0;
  cerr<<"this is step one "<<endl;
  i=0;
  cout<<"the a2 is "<<a2<<endl;
  while(a2[i] != '\0')
     i++;
  j=i;
  cerr<<"this is step two "<<endl;
  while(a2[j] != '/')
     j--;
  int z;
  j++;
  for(z=0;z<=i-j;z++)
    logfilename[z] = a2[z+j];
  logfilename[z] = '\0';
  strcat(logfilename,".log");
  cerr<<"the logfilename is "<<logfilename<<endl;
  strcpy(logstring,"date > ");
  strcat(logstring,logfilename);
  system(logstring);
  while(a3[argvint]!='\0')
  {
     argvint++;
  }
  if(a3[argvint-1]!='/')
  {
      strcat(a3,"/");
  }
  cout<<"a3 is "<<a3<<endl;
  sprintf(command, "mkdir  -p %s  ", a3 );
  system(command);
  strcpy(qfilename,a1);
  ifstream inputfile(qfilename,ios::in);
  if(!inputfile)
  {
      cout<<"the "<<qfilename<<" can't open "<<endl;
      exit(0);
  }
  inputfile.getline(quline,1000,'\n');
  checkNumberLine(quline); //Why is the query matrix limited to 10x10? *rds 11/08/2010
  inputfile.getline(quline,1000,'\n');
  qrow = quNumOfele((char *)quline , qElecol);
  cout<<"qrow is "<<qrow<<endl;
  quryMatrix = new char *[qrow];
  for(i=0;i<qrow;i++)
  {
     quryMatrix[i] = new char[qrow];
  }
  //this step I will get the quryMatrix.
    formqMatrix(quryMatrix, qrow,inputfile,totalhandqury,totalrequire);
  /*
  for(i=0;i<qrow;i++)
  {
    for(j=0;j<qrow;j++)
     cout<<quryMatrix[i][j];
    cout<<endl;
  }
*/
 //above step the qury matrix has been properly formed.
 //the follow step is to form the intaction matrix.
   strcpy(dirname , a2);
   ifstream inputfile1(dirname,ios::in);
   if(!inputfile1)
   {
       cout<<"the "<<dirname<<" file can't open "<<endl;
       exit(0);
   }
   readLincon = 0;
   countpid = 0;
   while(inputfile1.getline(intMline,10000000,'\n'))
   {
       //check = intMline;
       //cout<<"the check size is "<<check.length()<<endl;
       if(check.length() > 10000000)
       {
           cout<<"the length is wrong "<<endl;
           exit(0);
       }
       //cout<<"the check is "<<check<<endl;
       if(intMline[0] !='s')
       {
          bool isHeader = (strstr(intMline, ".ssd") != NULL);
          if(isHeader)
          {
               // A header (".ssd" line) begins a record. Discard any incomplete
               // pending record (header with no matrix) and its stray sheets, so
               // one malformed/orphan block can't desync the rest of the scan.
               if(judge1)
               {
                  for(i=0;i<mrow;i++)
                      delete [] interActionM[i];
                  delete [] interActionM;
               }
               totalsheet.clear();
               judge1 = false;
               judge2 = false;
               IntMnumEl.clear();
               intElecol.clear();
               mrow = intMnumofele(IntMnumEl,intMline,intElecol,pid ) ;
               if (mrow < 0) {
                   // Malformed DB block (orphan matrix or empty line where header
                   // expected). Reset state and try to recover sync by treating
                   // the next line as a fresh start; don't increment readLincon
                   // here so the next non-sheet line will retry as header.
                   judge1 = false;
                   judge2 = false;
                   continue;
               }
               cout<<"the step 1 "<<endl;
               cout<<"the mrow is "<<mrow<<endl;
              // cout<<"the intMline is "<<intMline<<endl;
               /*
               cout<<"the size of intElecol is "<<intElecol.size()<<endl;
               for(int w=0;w<intElecol.size();w++)
                  cout<<intElecol[w].geteleTyp();
               exit(0);
               */
               interActionM = new char *[mrow];
               for(i=0;i<mrow;i++)
               {
                  interActionM[i] = new char[mrow];
               }
               judge1 = true;
          }
          else if(judge1 && !judge2)
          //non-header, non-sheet line = the interaction matrix. Only valid when
          //a header is pending; an orphan matrix (no header) is skipped, which
          //keeps the scan in sync instead of silently dropping the rest.
          {
              getInterActionM(interActionM, mrow , intMline,totalsheet);
              judge2 = true;
          }
          if(judge1 == true && judge2 == true)
          {
		  cout<<"inside j1j2 loop"<<endl;
              //this step when the interAction matrix form finished, I just do the search.
              totalpass.clear();
              //this function I will print the strand with no h-bond with other strands
              //in the ../sheetbug/total.txt file
              checksheetH(interActionM , mrow , totalsheet,pid,intElecol);
              searchM(intElecol,interActionM, mrow,qElecol,quryMatrix,qrow,totalpass,IntMnumEl);
              //IntMnumEl this parameter is the matrixElment data type vector
              //totalhandqury,totalrequire these two parameter are the handness information and
              //other require information in the qury matrix.
              //totalsheet is the
	      cout<<"selectMatrix entry point"<<endl;
	      cout<<"selectMatrix: totalrequire size: "<<totalrequire.size()<<endl;
              selectMatrix(IntMnumEl,totalhandqury,totalrequire,totalpass,totalsheet,intElecol,interActionM);
	      cout<<"printOuptfile entry point"<<endl;
              printOuptfile(totalpass,IntMnumEl,pid,a3);
	      cout<<"the step 3"<<endl;
	      cout<<"totalsheet size: "<<totalsheet.size()<<endl;
              totalsheet.clear();
	      cout<<"the step 3.5"<<endl;
              for(i=0;i<mrow;i++)
                  delete [] interActionM[i];
              delete [] interActionM;
              judge1 = false;
              judge2 = false;
          }
           readLincon++;
       }
       else
       {
             //here I deal with the sheet part, I mean I read the sheet and form the sheet 
             //data structure
             int tmpint,bound;
             int index;
             int starposition = 18;
             temp = intMline;
             tempstring.assign(temp,6,6);
             tmpint = atoi(tempstring.c_str());
             if(tmpint != 0)
             {
                sheetptr = new sheet();
                tempstring.assign(temp,12,6);
                tmpint = atoi(tempstring.c_str());
                bound = tmpint;
                for(index=0;index<bound;index++)
                {
                   tempstring.assign(temp,starposition,6);
                   tmpint = atoi(tempstring.c_str());
                   //cout<<"the tmpint is "<<tmpint<<endl;
                   tmpint = tmpint - 1;
                   sheetptr->addelement(tmpint);
                   starposition = starposition + 6;
                }
                totalsheet.push_back(*sheetptr);
                delete sheetptr;
             }
        }
    }
  strcpy(logstring,"date >> ");
  strcat(logstring,logfilename);
  system(logstring);
  // Removed a per-query `grep ".ssd" <DB> > junk; wc -l junk` sanity count: it
  // re-read the entire 2.6GB DB once per query (6336 queries -> ~16TB of extra
  // reads across a full sweep) just to log a record count the program already
  // tracks as countpid.
  logptr = fopen(logfilename,"a");
  if(logptr == NULL)
  {
     cout<<"the file "<<logfilename<<" can't open "<<endl;
     exit(0);
  }
  fprintf(logptr,"pdb number counted in program: %d\n", countpid);
  fclose(logptr);
  cout<<"the countpid is "<<countpid<<endl;
}

bool searchControl::formajob(char *filelist,char *sendarray1,int &size)
{
    static int i=0;
    int j = 0;
    size = 0;
    int starnum = 0;
    while(filelist[i] !='\0')
    {
      sendarray1[j] = filelist[i];
      size++;
      i++;
      j++;
      if(filelist[i] == '*')
      {
          sendarray1[j] = '\0';
          i++;
          break;
      }
    }
    sendarray1[size] = '\0';
    if(size>0)
       return true;
    else
       return false;
}
      
//this filelistname should be replaced by the a directory
//thus function is used to check if the first line of the qury matrix is started with 1,2,3,4,5,6
void searchControl::checkNumberLine(char *qufile)
{
   int i = 0;
   if(qufile[0]!=' '&&qufile[0]!='1')
   {
      cout<<"in your qury file, the first line must be number index line"<<endl;
      cout<<"your number index line is : "<<qufile<<endl;
      exit(0);
   }
   if(qufile[0]==' ')
   {
      cout<<"in your qury file, the number index line, the first position can't be a space"<<endl;
      cout<<"your number index line is : "<<qufile<<endl;
      exit(0);
   }
   while(qufile[i+1]!='\0')
   {
      if(qufile[i]!=' '&&qufile[i+1]!=' ' && i < 10)
      {
          cout<<"your number index line, each number must be separated by a space"<<endl;
          cout<<"your number index line is : "<<qufile<<endl;
          //exit(0);
      }
      i++;
   }
}
    
bool searchControl::sameOrdiffsheet(vector<sheet> &totalsheet,vector<int>eleid)
{
    int i;
    int j;
    int u;
    bool same = true;
   /*
    for(i=0;i<eleid.size();i++)
      cout<<eleid[i]<<" ";
    cout<<endl;
   
    for(i=0;i<totalsheet.size();i++)
    {
        for(j=0;j<totalsheet[i].getSheet().size();j++)
        {
           cout<<totalsheet[i].getSheet()[j]<<" ";
        }
        cout<<endl;
    }
   */
    for(i=0;i<totalsheet.size();i++)
    {
       same = false;
       for(u=0;u<eleid.size();u++)
       {
          same = false;
          for(j=0;j<totalsheet[i].getSheet().size();j++)
          {
             if(eleid[u] == totalsheet[i].getSheet()[j])
             {
                 same = true;
             }   
          }
          if(same == false)
             break;
        }
      //cout<<"in loop same is "<<same<<endl; 
     if(same == true)
       break;
    }
    //cout<<"the same is "<<same<<endl;
    return same;
}

void searchControl::printgraph(vector<elecol> intEle ,int sou , int des ,vector<int> &path)
{
    int i;
    if(des == sou)
      path.push_back(sou);
    else if(intEle[des].getparent() == -1)
    { 
      //    cout<<"no path from "<<sou<<" to "<<des<<endl;
     //   exit(0);
          ;
    }
    else
    {
        printgraph( intEle , sou , intEle[des].getparent() ,path) ;
        path.push_back(des);
    }
}
      
int searchControl::paraOrantisearch(vector<elecol> &intEle,int source,int destination,vector<int>&colum,vector<int> ssheet,
                                    char **intM)
{
    int i ,j , z;
    int soucolum , descolum;
    int soucolum1 , descolum1;
    vector<elecol> onelepass;
    vector<int> neighbor;
    queue<elecol> que;
    elecol head;
    vector<int> path;
    bool judge = false;
    int headindex = 0;
    if(source >=colum.size() || destination >= colum.size())
    {
        cout<<"the in paraOrantisearch function source or destination or colum size data is wrong"<<endl;
        exit(0);
    }
    /*
    cout<<"the colum is "<<endl;
    for(i=0;i<colum.size();i++)
      cout<<colum[i]<<" ";
    cout<<endl;
    cout<<"the ssheet content is "<<endl;
    for(i=0;i<ssheet.size();i++)
      cerr<<ssheet[i]<<" ";
    cout<<endl;

    cout<<"the intEle size is "<<intEle.size()<<endl;
    for(i=0;i<intEle.size();i++)
    {
       cout<<"No "<<i<<" content is "<<endl;
       intEle[i].print();
    }
    */
    for(i=0;i<ssheet.size();i++)
    {
       onelepass.push_back(intEle[colum[ssheet[i]-1]]);
    }
    
    for(i=0;i<onelepass.size();i++)
    {
        neighbor = onelepass[i].getneighbor();
        for(j=0;j<neighbor.size();j++)
        {
            judge = false;
            for(z=0;z<onelepass.size();z++)
            {
		if(neighbor[j] == onelepass[z].getcol())
                   judge = true;
            }
            if(judge == false)
               onelepass.push_back(intEle[neighbor[j]]);
        }
        neighbor.clear();
    } 
    bool findparent = false;
    soucolum = colum[source];
    descolum = colum[destination];
    for(i=0;i<onelepass.size();i++)
    {
        if(onelepass[i].getcol() == soucolum)
           soucolum1 = i;
        if(onelepass[i].getcol() == descolum)
           descolum1 = i;
    }
    soucolum = soucolum1;
    descolum = descolum1;
    for(i=0 ;i<onelepass.size();i++)
      onelepass[i].setcolorparent("white" , -1);
    //cout<<"the soucolum is "<<soucolum<<endl;
    //cout<<"the descolum is "<<descolum<<endl;
    onelepass[soucolum].setcolorparent("gray" , -1);
    que.push(onelepass[soucolum]);
    bool findnei = false;
    while(que.size() != 0)
    {
       head = que.front();
       findparent = false;
       for(i=0;i<onelepass.size();i++)
       {
           if(head.getcol() == onelepass[i].getcol())
           {
               headindex = i;
               findparent = true;
           }
       }
       if(findparent == false)
       {
          cout<<"in garaph search something is wrong"<<endl;
          exit(0);
       }
       neighbor = head.getneighbor();
       for(j=0;j<neighbor.size();j++)
       {
          findnei = false;
          for(i=0;i<onelepass.size();i++)
          {
     //        cout<<"the neighbor[j] is "<<neighbor[j]<<endl;
             if(neighbor[j] == onelepass[i].getcol())
             {
                findnei = true;
                if(onelepass[i].getcolor() == "white")
                {
                   onelepass[i].setcolorparent("gray" ,headindex);
                   que.push(onelepass[i]);
                }
             }
          }
          if(findnei == false)
          {
             cout<<"the neighbor data or onelepass data wrong "<<endl;
             exit(0);
          }
       }
       que.pop();
    } 
    path.clear();
    printgraph(onelepass ,soucolum , descolum ,path);

    vector<char> path1;
    int row1 , col1;
    char calculate;
    int judgep = 2;
    //0 means parallel , 1 means antiparall ,2 means they aren't parallel or antiparallel;
    if(path.size() == 1)
    {
        cout<<"something is wrong in judge parallel or antiparallel function"<<endl;
        exit(0);
    }
    else
    {
        if(path.size()>1)
        {
           for(i=0;i<path.size()-1;i++)
           {
              row1 = onelepass[path[i]].getcol();
              col1 = onelepass[path[i+1]].getcol();
              path1.push_back(intM[row1][col1]);
           }
        }
          
    }
    /*
    cout<<"path1 size is "<<path1.size()<<endl;
    for(i=0 ;i<path1.size();i++)
      cout<<path1[i]<<" ";
    cout<<endl;
    */
    if(path1.size() == 0)
       judgep = 2;
    if(path1.size() == 1)
    {
       if(path1[0] == 'c')
          judgep = 0; 
       else if(path1[0] == 't')
          judgep = 1; 
       else
       {
          cout<<"in paraOr function path1 data is wrong "<<endl;
          exit(0);
       }
    }
    if(path1.size()>1)
    {
        calculate = path1[0];
        for(i=1;i<path1.size();i++)
        {
            if(calculate == 'c' && path1[i] == 'c')
               calculate = 'c';
            if(calculate == 'c' && path1[i] == 't')
               calculate = 't';
            if(calculate == 't'&& path1[i] == 'c')
               calculate = 't';
            if(calculate == 't' && path1[i] == 't')
               calculate == 'c';
        }
        if(calculate == 'c')
          judgep = 0;
        else if(calculate == 't')
          judgep = 1;
        else
        {
          cout<<"in paraOr function path1size great than 1 data is wrong "<<endl;
          exit(0);
        }
    } 
    
    return judgep;
}
//the parameter b is the vector of data type matrixElment , oneline is the string that got from
//the yioutputmatrix, intEle is the vector of the data type elecol, this data type is uesd
//for the later matrix search.
void searchControl::selectMatrix(vector<matrixElment> &a,vector<handness> &hand,vector<require> &tR,vector<fpass>&allps,vector<sheet> &totalsheet, vector<elecol> &intEle,char **intM)
{
   fpass temppass;
   vector<fpass> last;
   int i=0;
   int j =0;
   int z =0;
   int u = 0;
   bool handbool;
   bool requirebool;
   int passIndex =0;
   vector<matrixElment> onematrix;
   matrixElment seg1;
   matrixElment seg2;
   matrixElment seg3;
   vector<int> colum;
   require tmpr;
   vector<int> QsheeteleId; 
   char handchar;
   int h0,h1,h2;
   int tempcolum; 
   int sheetint =0;
   vector<int> ssheet;
   //first step I will get the matrix at first.
   cout<<"this is selectMatrix "<<endl;
   cout<<"the tR size is "<<tR.size();
   //cout<<"the tmpr.getrequireType() is "<<tR[u].getrequireType()<<endl;
   cout<<"the totalpasee size is "<<allps.size()<<endl;
   
   for(i=0;i<allps.size();i++)
     allps[i].print();
   cout<<"the allps size is "<<allps.size()<<endl;
   for(i=0;i<allps.size();i++)
   {
      cout<<"this is No +++++++++++++++ "<<i<<endl;
      allps[i].print();
   }
   //exit(0);
   
   for(passIndex=allps.size()-1;passIndex>=0;passIndex--)
   {
       temppass = allps[passIndex];
       if(tR.size()==0 && hand.size()==0 )
          last.push_back(temppass);
       colum.clear();
       for(j=0;j<temppass.getparent().size();j++) 
       {
         tempcolum = temppass.getparent()[j].getcol();
         colum.push_back(tempcolum); 
       }
       onematrix.clear();
       for(u=0;u<colum.size();u++)
       {
          onematrix.push_back(a[colum[u]]);
       }
       cout<<"the chain is "<<onematrix[3].getChainName();
       handbool = true;
       for(i=0;i<hand.size();i++)
       {
           handchar = hand[i].gethandness();
           if(hand[i].getIdgroup().size()<3)
           {
              cout<<"the qury file handness qury is wrong "<<endl;
              exit(0);
           }
           char rechar;
           h0 = hand[i].getIdgroup()[0];
           h1 = hand[i].getIdgroup()[1];
           h2 = hand[i].getIdgroup()[2];
           seg1 =  onematrix[h0-1];
           seg2 =  onematrix[h1-1]; 
           seg3 =  onematrix[h2-1];
           if(h0>onematrix.size()||h1>onematrix.size()||h2>onematrix.size())
           {
               cout<<"in qury matrix h0 or h1 or h2 is wrong "<<endl;
               exit(0);
           }
           rechar = chirality(seg1,seg2,seg3);
           cout<<"the rechar is "<<rechar<<endl;
           cout<<"the handchar is "<<handchar<<endl;
           if(rechar != handchar)
           {
              handbool = false;
              break;
           }
       }
       //this step I will judge the lengh require;
       if(handbool == false)
          requirebool = false;
       else
          requirebool = true;
       
       int eleid,length;
       if(handbool == true)
       {
          for(u=0;u<tR.size();u++)
          {
             tmpr = tR[u]; 
             cout<<"the tmpr.getrequireType() is "<<tR[u].getrequireType()<<endl;
             requirebool = true;
             if(tmpr.getrequireType() == "length")
             {
                eleid = tmpr.geteleId(); //need change
                length = onematrix[eleid - 1].getlenght();//need change
                cout<<"the length is "<<length<<endl;
                cout<<"tmpr.getstartRange() "<<tmpr.getstartRange()<<endl;
                cout<<"tmpr.getendRange()   "<<tmpr.getendRange()<<endl;
                if(length<tmpr.getstartRange() || length>tmpr.getendRange())
                {
                   requirebool = false;
                   break;
                }
             }
             int sheetid;
             if(tmpr.getrequireType() == "sheetD")
             {
                QsheeteleId.clear();
                for(sheetint=0;sheetint<tmpr.getElementSet().size();sheetint++)
                {
                   QsheeteleId.push_back(colum[tmpr.getElementSet()[sheetint]-1]);
                }
                if(sameOrdiffsheet(totalsheet,QsheeteleId)==true)
                {
                   requirebool = false;
                   break;
                } 
                else 
                {
                   requirebool = true;
                }
             }
             if(tmpr.getrequireType() == "sheetS")
             {
                QsheeteleId.clear();
                ssheet.clear();
                ssheet = tmpr.getElementSet();
                for(sheetint=0;sheetint<tmpr.getElementSet().size();sheetint++)
                {
                   QsheeteleId.push_back(colum[tmpr.getElementSet()[sheetint]-1]);
                }
               
                for(int v=0;v<QsheeteleId.size();v++)
                {
                    cout<<"the sheet index is "<<QsheeteleId[v]<<"  ";
                }
                cout<<endl;
                
                if(sameOrdiffsheet(totalsheet,QsheeteleId)==true) 
                {
                    requirebool = true;
                }
                else
                {
                    requirebool = false;
                    break;
                }
             }
             if(tmpr.getrequireType() == "chainAll")
                requirebool = true;
             if(tmpr.getrequireType() == "chainS")
             {
                int u=0;
                int elementId;
                char cname;
                if(tmpr.getElementSet().size()<2)
                {
                    cout<<"qury file chainS line is wrong "<<endl;
                    exit(0);
                }
                elementId = tmpr.getElementSet()[0]-1;
                cname = onematrix[elementId].getChainName();   
                cout<<"the cname is "<<cname<<endl;
                cout<<"the mpr.getElementSet().size() "<<tmpr.getElementSet().size()<<endl;
                for(u=0;u<tmpr.getElementSet().size() ;u++)
                {
                    elementId = tmpr.getElementSet()[u]-1; 
                    cout<<"the elementId is "<<elementId<<endl;
                    if(cname!=onematrix[elementId].getChainName())
                    {
                        requirebool = false;
                        break;
                    }
                    else
                    {
                    	cout<<"it is eaqula "<<endl;
                        requirebool = true;
                    }
                 }
                 if(requirebool == false)
                    break;
                 
             }    
              
             if(tmpr.getrequireType() == "chainD")
             {
                int u=0;
                int elementId;
                char cname;
                if(tmpr.getElementSet().size()<2)
                {
                    cout<<"qury file chainS line is wrong "<<endl;
                    exit(0);
                }
                elementId = tmpr.getElementSet()[0]-1;
                cname = onematrix[elementId].getChainName();   
                cout<<"the cname is "<<cname<<endl;
                cout<<"the mpr.getElementSet().size() "<<tmpr.getElementSet().size()<<endl;
                for(u=0;u<tmpr.getElementSet().size() ;u++)
                {
                    elementId = tmpr.getElementSet()[u]-1; 
                    cout<<"the elementId is "<<elementId<<endl;
                    cout<<"onematrix[elementId].getChainName() "<<onematrix[elementId].getChainName()<<endl;
                    if(cname!=onematrix[elementId].getChainName())
                    {
                        requirebool = true;
                        break;
                    }
                    else
                    {
                        requirebool = false;
                    }
                }
                 //cout<<"the requirebool in chain section "<<requirebool<<endl;
                 if(requirebool == false)
                    break;
                 
             }    
             if(tmpr.getrequireType() == "parallel" || tmpr.getrequireType() == "antiparallel")
             {
                int u=0;
                int t = 0;
                int elementId1,elementId2;
                char cname;
                vector<int> para;
                point aF , aB , bF , bB;
                double angle1;
                int judg;
                
                for(t=0 ;t<intEle.size();t++)
                   intEle[t].print();
                cout<<endl;
                for(z=0;z<colum.size();z++)
                   cout<<colum[z]<<" ";
                cout<<endl;
                exit(0);
                
                if(tmpr.getElementSet().size()!=2)
                {
                    cout<<"qury file parallel is only 2 choice "<<endl;
                    exit(0);
                } 
                for(u=0;u<tmpr.getElementSet().size() ;u++)
                {
                    //here I will call the graph search,then use the element id as source and destination
                    //then get the path.
                    if(u==0)
                    {
                       elementId1 = tmpr.getElementSet()[u]-1;
                       aF = onematrix[elementId1].getFpoint();
                       aB = onematrix[elementId1].getLpoint();
                    }
                    else if(u==1)
                    {
                       elementId2 = tmpr.getElementSet()[u]-1;
                       bF = onematrix[elementId2].getFpoint();
                       bB = onematrix[elementId2].getLpoint();
                    }
                    else
                    {
                       cout<<"in parallel u value is wrong "<<endl;
                       exit(0);
                    }
                }
                 judg = paraOrantisearch(intEle,elementId1,elementId2,colum,ssheet,intM);
                 if(tmpr.getrequireType() == "parallel")
                 {
                    if(judg !=0)
                    { 
                        requirebool = false;
                        break;
                    }
                    else
                        requirebool = true;
                 }
 
                 if(tmpr.getrequireType() == "antiparallel")
                 {
                    if(judg != 1)
                    {
                       requirebool = false;
                       break;
                    }
                    else
                       requirebool = true;
                 }
             }
         }
         if(requirebool == true && hand.size() !=0 || requirebool == true && tR.size() !=0)
         {
              last.push_back(temppass);
         }
          
      }
   }   
  
   allps.clear();
   for(int p=0;p<last.size();p++)
     allps.push_back(last[p]);
   cout<<"the allps is "<<allps.size()<<endl;
}
//this function is used to judge the charility of the three elements
char searchControl::chirality(matrixElment seg1,matrixElment seg2,matrixElment seg3)
{
    double N1x, N1y, N1z;
    double N2x, N2y, N2z;
    double N3x, N3y, N3z;
    double Nx, Ny, Nz;
    double Mx, My, Mz;
    double ang;
    double m1x,m1y,m1z;
    double m2x,m2y,m2z;
    double m3x,m3y,m3z;
    point seg1initPoint,seg1endPoint,seg2initPoint,seg2endPoint;
    point seg3endPoint,seg3initPoint;
    seg1initPoint = seg1.getFpoint();
    seg1endPoint  = seg1.getLpoint();
    seg2initPoint = seg2.getFpoint();
    seg2endPoint  = seg2.getLpoint();
    seg3initPoint = seg3.getFpoint();
    seg3endPoint  = seg3.getLpoint();
    N1x = seg1endPoint.xcoord-seg1initPoint.xcoord;
    N1y = seg1endPoint.ycoord-seg1initPoint.ycoord;
    N1z = seg1endPoint.zcoord-seg1initPoint.zcoord;
    ang = angle(seg1initPoint.xcoord,seg1initPoint.ycoord,seg1initPoint.zcoord,
               seg1endPoint.xcoord,seg1endPoint.ycoord,seg1endPoint.zcoord,
               seg2initPoint.xcoord,seg2initPoint.ycoord,seg2initPoint.zcoord ,
               seg2endPoint.xcoord,seg2endPoint.ycoord,seg2endPoint.zcoord);
    if (ang<90.001 )
        {
                N2x = seg2endPoint.xcoord-seg2initPoint.xcoord;
                N2y = seg2endPoint.ycoord-seg2initPoint.ycoord;
                N2z = seg2endPoint.zcoord-seg2initPoint.zcoord;

        }
        else
        {
                N2x = seg2initPoint.xcoord - seg2endPoint.xcoord;
                N2y = seg2initPoint.ycoord - seg2endPoint.ycoord;
                N2z = seg2initPoint.zcoord - seg2endPoint.zcoord;
        }
        ang = angle(seg1initPoint.xcoord,seg1initPoint.ycoord,seg1initPoint.zcoord,
                    seg1endPoint.xcoord,seg1endPoint.ycoord,seg1endPoint.zcoord,
                    seg3initPoint.xcoord,seg3initPoint.ycoord,seg3initPoint.zcoord,
                    seg3endPoint.xcoord,seg3endPoint.ycoord,seg3endPoint.zcoord);
        if (ang< 90.001)
        {
                N3x = seg3endPoint.xcoord-seg3initPoint.xcoord;
                N3y = seg3endPoint.ycoord-seg3initPoint.ycoord;
                N3z = seg3endPoint.zcoord-seg3initPoint.zcoord;
        }
        else
        {
                N3x = seg3initPoint.xcoord - seg3endPoint.xcoord;
                N3y = seg3initPoint.ycoord - seg3endPoint.ycoord;
                N3z = seg3initPoint.zcoord - seg3endPoint.zcoord;
        }

        Nx = ( N1x + N2x + N3x ) / 3;
        Ny = ( N1y + N2y + N3y ) / 3;
        Nz = ( N1z + N2z + N3z ) / 3;

#ifdef HelixSearch
        Nx = seg1endPoint.xcoord-seg1initPoint.xcoord; 
/*use the direction of the first element for Helix Search */
#endif
    m1x = (seg1endPoint.xcoord+seg1initPoint.xcoord)/2;
    m1y = (seg1endPoint.ycoord+seg1initPoint.ycoord)/2;
    m1z = (seg1endPoint.zcoord+seg1initPoint.zcoord)/2;
    m2x = (seg2endPoint.xcoord+seg2initPoint.xcoord)/2;
    m2y = (seg2endPoint.ycoord+seg2initPoint.ycoord)/2;
    m2z = (seg2endPoint.zcoord+seg2initPoint.zcoord)/2;
    m3x = (seg3endPoint.xcoord+seg3initPoint.xcoord)/2;
    m3y = (seg3endPoint.ycoord+seg3initPoint.ycoord)/2;
    m3z = (seg3endPoint.zcoord+seg3initPoint.zcoord)/2;
    Mx = (m2y-m1y)*(m3z-m1z)-(m2z-m1z)*(m3y-m1y);
    My = (m2z-m1z)*(m3x-m1x)-(m2x-m1x)*(m3z-m1z);
    Mz = (m2x-m1x)*(m3y-m1y)-(m2y-m1y)*(m3x-m1x);
    ang = angle(0,0,0,Nx,Ny,Nz,0,0,0,Mx,My,Mz);
    fprintf(stdout,"angle in chirality: %f %5s %5s %5s\n",ang,(seg1.getBeginNum()).c_str(),(seg2.getBeginNum()).c_str(),(seg3.getBeginNum()).c_str());
    if(ang<90){
       return 'R';
       //cout<<ang<<" "<<seg1.getBeginNum()<<" "<<seg2.getBeginNum()<<" "<<seg3.getBeginNum();
    }  
    else if(ang>90.001){
       return 'L';
       //cout<<ang<<" "<<seg1.getBeginNum()<<" "<<seg2.getBeginNum()<<" "<<seg3.getBeginNum();
    }   
    else{
    	 //cout<<ang<<" "<<seg1.getBeginNum()<<" "<<seg2.getBeginNum()<<" "<<seg3.getBeginNum();
       return ' ';/* change 'N' to ' '*/
    }   
}

void searchControl::getInterActionM(char ** intM,int mrow , char *oneline,vector<sheet> &totalsheet)
{
    int i,j;
    int index = 0;
    sheet *shptr;
    string lengthstr;
    cout<<"getInterActionM entry point "<<endl;
    cout<<"the mrow is "<<mrow<<endl;
    lengthstr = oneline;
    cout<<"the length oneline is "<<lengthstr.length()<<endl;
    cout<<"the mrow is "<<mrow<<endl;
    for(i=0;i<mrow;i++)
       for(j=0;j<mrow;j++)
       {
         intM[i][j] = ' ';
       }
    cout<<"this is output 3.5"<<endl;
    i = -1;
    j = 0;
    while(oneline[index] != '\0')
    {
       if(oneline[index] == '*')
       {
           i++;
           // Guard against matrix-line longer than expected (more '*' markers than mrow).
           // Without this, intM[i][j] writes past the heap allocation of mrow rows.
           if (i >= mrow) break;
           j = i;
           intM[i][j] = oneline[index];
           j++;
       }
       else
       {
           // Same guard for row-internal characters that would exceed the row width.
           if (i < 0 || i >= mrow || j >= mrow) { index++; continue; }
           intM[i][j] = oneline[index];
           j++;
       }
       index++;
    }

    for(i=0;i<mrow;i++)
    {
        for(j=0;j<mrow;j++)
          cout<<intM[i][j]<<" ";
        cout<<endl;
    }
 
    for(i=0;i<mrow;i++)
       for(j=i;j<mrow;j++)
    {
       intM[j][i] = intM[i][j];
    }

    for(i=0;i<mrow;i++)
    {
       for(j=0;j<mrow;j++)
          cout<<intM[i][j];
       cout<<endl;
    }

    cout<<"out of function "<<endl;
}
//this function is to initialize the varible b and intEle, 
//matricElment this data structure contain the element information,come from the generate matrix
//element type and range and start and end point and their coordinates

int searchControl::intMnumofele(vector<matrixElment> &b  , char *oneline,vector<elecol> &intEle,char *pid )
{
    int i , j;
    int numberElment;
    string beRnum, endRnum;
    int thiseleLength;
    int size;
    elecol *eleptr;
    point a;
    char chainName, eletype;
    double xcoord,ycoord,zcoord;
    string cpline;
    matrixElment *mptr;
    string temp;
    int beginposition = 36;
    cpline = oneline;
    size = cpline.size();
    //this maybe a potential bug "
    if(cpline.size()<14)
    {
       cout<<"the yimatrixoutput file is wrong, please check "<<endl;
       exit(0);
    }
    // Validate that this looks like a header line. The discriminator is the
    // ".ssd" filename suffix — header lines have it (entry IDs end in .ssd),
    // orphan matrix lines do not (they are packed interaction codes like
    // "*u---v---*..."). The original strict digit-prefix check rejected AFDB
    // headers (named e.g. "dpam_A0A011QYY6_nD2.ssd") and made AFDB DB scans
    // return 0 hits immediately. Just checking for .ssd works for both
    // ECOD-style numeric IDs and AFDB-style alphanumeric IDs.
    if (cpline.find(".ssd") == string::npos)
    {
       cerr << "intMnumofele: non-header line detected ("
            << cpline.substr(0, 40) << "...) -- DB has malformed block; skipping" << endl;
       return -1;
    }
    // pid is a 20-byte buffer in oneprocess() (line 195). Copying 31 chars overflows
    // it. Read at most 19 chars (leaving room for the null terminator).
    temp.assign(cpline,0,19);
    strncpy(pid, temp.c_str(), 19);
    pid[19] = '\0';
    countpid++;
    int u=0;
    while(u < 19 && pid[u]!=' ')
      u++;
    pid[u] = '\0';
    cout<<"the pid is "<<pid<<endl;
    temp.assign(cpline,32,4);
    numberElment = atoi(temp.c_str());
    int col = 0;
    while(beginposition<size)
    {
        mptr = new matrixElment();
        temp.assign(cpline,beginposition,1);
        beginposition = beginposition + 1;
        eletype = temp[0];
        temp.assign(cpline,beginposition,1);
        chainName = temp[0];
        beginposition = beginposition + 1;
        beRnum=temp.assign(cpline,beginposition,5);//need my change
        //beRnum = atoi(temp.c_str());
        //cout<<"the beRnum is "<<beRnum<<endl;
        beginposition = beginposition + 7;
        endRnum=temp.assign(cpline,beginposition,5);
        //endRnum = atoi(temp.c_str());//need my change
        //cout<<"the endRnum is "<<endRnum<<endl;
        beginposition = beginposition + 6;
        temp.assign(cpline,beginposition,4);
        thiseleLength=atoi(temp.c_str());
        beginposition=beginposition+5;
        eleptr = new elecol(eletype , col,chainName);
        intEle.push_back(*eleptr);
        delete eleptr;
        mptr->setBeginNum(beRnum);
        mptr->setEndNum(endRnum);
        mptr->setElmentType(eletype);
        mptr->setChainName(chainName);
        mptr->setlength(thiseleLength);
        for(i=0;i<2;i++)
        {
           temp.assign(cpline,beginposition,8);
           beginposition = beginposition + 8;
           xcoord = atof(temp.c_str());
           temp.assign(cpline , beginposition ,8);
           ycoord = atof(temp.c_str());
           beginposition = beginposition + 8;
           temp.assign(cpline,beginposition,8);
           zcoord = atof(temp.c_str());
           a.xcoord = xcoord; a.ycoord = ycoord ; a.zcoord = zcoord;
           /*
           cout<<"the xcoord is "<<xcoord<<endl;
           cout<<"the ycoord is "<<ycoord<<endl;
           cout<<"the zcoord is "<<zcoord<<endl;
           */
           if(i==0)
           {
             mptr->setFpoint(a);
           } 
           if(i==1)
           {
             mptr->setLpoint(a);
           }
           beginposition = beginposition + 8;
        }
        b.push_back(*mptr);
        delete mptr;
        col++;
      }
    /*
     for(i=0;i<b.size();i++)
        b[i].print();
     for(i=0;i<intEle.size();i++)
        intEle[i].print();
    */
     return b.size();
}     
//this function I read the data from the qury matrix file and get the
//qury matrix.and also get the handness character in the file.

void searchControl::formqMatrix(char ** matrix,int row,ifstream &a,vector<handness> &hand ,vector<require> &tR) 
{
     int i , j;
     handness temphand;
     int g=0;
     int w = 0;
     int num = 0;
     int u =0;
     int lineNum = 0;
     int firstEl,secondEl,thirdEl;
     char handChar;
     char line[1000];
     char temchar[10];
     require *rptr;
     int   eleid,start,end;
     string eletype;
     string temp;
     string wholeLine;
     char *tokenPtr,*tokenPtr1;
     int bposition;
     char line1[500];
     bool firstime = false;
     int checkspace=0;
     for(i=0;i<row;i++)
        for(j=0;j<row;j++)
        matrix[i][j] = ' ';
     i=0;
     j=0;

     while(a.getline(line ,1000,'\n'))
     {
        lineNum++;
        if(lineNum<=row)
        {
           //cout<<"the line is "<<line<<endl;
           checkspace =0;
           while(line[checkspace+1]!='\0')
           {
               if(line[checkspace]!=' '&& line[checkspace+1]!=' ')
               {
                    cout<<"in your qury matrix each column must be separated by only one space"<<endl;
                    cout<<"please check your qury matrix"<<endl;
                    exit(0);
               }
               checkspace++;
           }
           w = 0;
           cout<<"line is "<<line<<endl;
           u = num;
           while(line[w]!='\0')
           {
               if(line[w] != ' ')
               {
                  matrix[i][u] = line[w];
                  u++;
               }
               w++;
            }
           i++;
           num++;
           if(num > row || i>row)
           {
               cout<<"in searchControl::formqMatrix function num wrong "<<endl;
               exit(0);
           }

        }
       if(lineNum > row )
       {
         wholeLine = line;
         strcpy(line1,line);
         cout<<"the wholeLine is "<<line<<endl;
         cout<<"the line1 is     "<<line1<<endl;
         if(wholeLine.find("handedness") == 0)
         {
            g=0;
            tokenPtr1 = strtok(line1," "); 
            while(tokenPtr1 !=NULL)
            {
                tokenPtr1 = strtok(NULL," ");
                g++;
            }
            cout<<"g value is "<<g<<endl;
            if(g<5 || g>5)
            {
                cout<<"g value is "<<g<<endl;
                cout<<"your handedness qury: "<<wholeLine<<" : is wrong"<<endl;
                exit(0);
            }

            g=0;
            tokenPtr = strtok(line," ");
            while(tokenPtr !=NULL)
            {
                tokenPtr = strtok(NULL," ");
                g++;
                if(g<=3)
                {
                    firstEl = atoi(tokenPtr);
                    if(firstEl > row)
                    {
                        cout<<"in your handedness qury index of your element is wrong"<<endl;
                        cout<<"your index of the element is "<<firstEl<<endl;
                        cout<<"it should be less or eaqual to "<<row<<endl;
                        exit(0);
                    }
                    temphand.addElmentId(firstEl);
                }
                if(g==4)
                {
                    temphand.sethandness(tokenPtr[0]);
                    cout<<"the tokenPtr[0] is "<<tokenPtr[0]<<endl;
                    hand.push_back(temphand);
                    temphand.clear();
                }
             }
          }

          else if(wholeLine.find("length") == 0)
          {
            g=0;
            tokenPtr1 = strtok(line1," ");
            rptr = new require();
            rptr->setRequireType("length");
            while(tokenPtr1 !=NULL)
            {
                tokenPtr1 = strtok(NULL," ");
                g++;
            }
           cout<<"g value is "<<g<<endl;
            if(g<5 || g>5)
            {
                cout<<"g value is "<<g<<endl;
                cout<<"your handedness qury: "<<wholeLine<<" : is wrong"<<endl;
                cout<<"it should be : length 1 E 2 4"<<endl;
                cout<<"1 is the index of element, E is the element type ,"<<endl;
                cout<<"2 4 means the length of this element is between 2 and 4, include 2 and 4"<<endl;
                exit(0);
            }
            g=0;
            tokenPtr = strtok(line," ");
            while(tokenPtr !=NULL)
            {
                tokenPtr = strtok(NULL," ");
                g++;
                if(g==1)
                {
                    eleid = atoi(tokenPtr);
                    if(eleid > row)
                    {
                        cout<<"in your length qury index of your element is wrong"<<endl;
                        cout<<"your index of the element is "<<eleid<<endl;
                        cout<<"it should be less or eaqual to "<<row<<endl;
                        exit(0);
                    }
                    rptr->setEleId(eleid);
                }
                if(g==2)
                {
                    rptr->setEleType(tokenPtr);
                }
                if(g==3)
                {
                    start = atoi(tokenPtr);
                }
                if(g==4)
                {
                    end = atoi(tokenPtr);
                    cout<<"the start is "<<start<<" the end is "<<end<<endl;
                    rptr->setStartEndRange(start, end);
                    tR.push_back(*rptr);
                    delete rptr;
                }
             }
       }
       else if(wholeLine.find("sheetD") == 0 ||wholeLine.find("sheetS") == 0)
       {
          rptr = new require();
          if(wholeLine.find("sheetD") == 0)
             rptr->setRequireType("sheetD");
          else
          rptr->setRequireType("sheetS");
          tokenPtr = strtok(line," ");
          while(tokenPtr !=NULL)
          {
              tokenPtr = strtok(NULL," ");
              if(tokenPtr !=NULL)
              {
                eleid = atoi(tokenPtr);
                if(eleid > row)
                {
                   cout<<"in your sheetD qury or sheetS qury the index of element is wrong"<<endl;
                   cout<<"your index of the element is "<<eleid<<endl;
                   cout<<"it should be less or eaqual to "<<row<<endl;
                   exit(0);
                }
                rptr->addElement(eleid);
              }
          }
         tR.push_back(*rptr);
         delete rptr;
      }
      else if(wholeLine.find("chainS")==0 || wholeLine.find("chainD")==0 
              || wholeLine.find("antiparallel") == 0 || wholeLine.find("parallel") == 0)
      {
          if(wholeLine.find("chainD")==0)
            chainjudge = false;
          rptr = new require();
          if(wholeLine.find("chainS")==0)
            rptr->setRequireType("chainS");
          if(wholeLine.find("chainD")==0) 
            rptr->setRequireType("chainD");
          if(wholeLine.find("parallel")==0)
            rptr->setRequireType("parallel");
          if(wholeLine.find("antiparallel")==0)
            rptr->setRequireType("antiparallel");
          tokenPtr = strtok(line," ");
          while(tokenPtr !=NULL)
          {
              tokenPtr = strtok(NULL," ");
              if(tokenPtr !=NULL)
              {
                eleid = atoi(tokenPtr);
                eleid = atoi(tokenPtr);
                if(eleid > row)
                {      
                   cout<<"in your chainD qury or chainS qury the index of element is wrong"<<endl;
                   cout<<"your index of the element is "<<eleid<<endl;
                   cout<<"it should be less or eaqual to "<<row<<endl;
                   exit(0);
                }
                rptr->addElement(eleid);
              }
         }
         tR.push_back(*rptr);
         delete rptr;
      }
      else if(wholeLine.find("chainAll") == 0)
      {
         chainjudge = false;
         rptr = new require();
         rptr->setRequireType("chainAll");
         tR.push_back(*rptr);
      }
      else{
          cout<<"the qury: "<<wholeLine<< ": is worng "<<endl;
          cout<<"we just have : handedness,sheetS,sheeD,length,chainS,chainD,length,parallel,antiparallel : these 7 types"<<endl;
          cout<<"make sure your spell is right"<<endl;
          exit(0);
      }
    }
 }
    /*
    for(i=0;i<hand.size();i++)
    {
        hand[i].print();
    }
    for(i=0;i<tR.size();i++)
    {
        tR[i].print();
    }
    */
    for(i=0;i<row;i++)
    {
        if(matrix[i][i]!='*')
        {
            cout<<"in your qury matrix, the "<<i+1<<" row,"<<i+1<<" colum is wrong "<<endl;
            cout<<"it should be character * "<<endl;
            cout<<"please check"<<endl;
            exit(0);
        }
    }
    char z;
    for(i=0;i<row;i++)
       for(j=0;j<row;j++)
       {
           z = matrix[i][j];
           if(z!='X'&& z!='x'&&z!='T'&&z!='C'&&z!='-'&&z!='t'&&z!='t'&&z!='v'&&z!='u'&&z!='N'&&z!='c'&&z!='*'&&z!=' ')
           {
               cout<<"character "<<z<<" is not in our defined character set"<<endl;
               cout<<"please check "<<i+1<<" row and "<<j+1<<"colum character"<<endl;
               exit(0);
           }
       }
}

//this function we get the content of the data structure of the quelecol.     
//the quelecol is store the information in the qury matrix, and store the element type and
//its correspond column
int searchControl::quNumOfele(char *line , vector<elecol> &quelecol)
{
   elecol *eptr;
   int i = 0;
   int j=0;
   int size;
   //I set the default chain is A;
   char chaiNa = 'A';
   if(line[i] == ' ')
   {
      cout<<"the first position in the element line must be H,E,X, it can't be space"<<endl;
      exit(0);
   }
   while(line[i+1]!='\0')
   {
      if(line[i]!=' '&&line[i+1]!=' ')
      {
         cout<<"each element type must be separated by one space"<<endl;
         cout<<"please check your qury matrix secondary element line,make sure each element is separated by space"<<endl;
         exit(0);
      }
      i++;
   }
   i=0;
   while(line[i]!='\0')
   {
      if(line[i] !='H'&&line[i] !='E'&&line[i] !='X'&&line[i] !=' '&&line[i] != 'L')
      {
         cout<<line[i]<<" is not in element types range "<<endl;
         cout<<"we just have four secondary element types: H,E,X,L"<<endl;
         cout<<"please check your qury matrix secondary element line,make sure you wrote the right character"<<endl;
         exit(0);
      }
      if(line[i] != ' ')
      {
          j++;
          eptr = new elecol(line[i] ,j,chaiNa);
          quelecol.push_back(*eptr);
          delete eptr;
      }
      i++;
   }
   size = quelecol.size();
   return size;
}
     
   
//vector<elecol> a is the parameter of the element of the secodary structure
//intM is the big interaction of the matrix, row1 is the number of the row
//because of this matrix is the square.
//!!!!!!!we have some different types of the requirements in the search algorithm
// (1) is the length reqirement for example: length 1 E 4 8,this means that element 1 in qury
//matrix the element type is E and length is greater than 4 and less than 8
// (2) is the sheet requirements for example: sheet same 1 2 3 4,means the elements 1 2 3 4 in qury 
//matrix should be in the same sheet.
//sheet different 1 4,means in qury matirx the element 1 and 4 are in the different sheet
void searchControl::searchM(vector<elecol>a,char ** intM, int row1,vector<elecol>b,char** quM,int row2,
                           vector<fpass> &totalpass,vector<matrixElment> &totalele)
{
    int i,j;
    //onepass is the structure of the que to keep all the possible ways of the combination
    queue<fpass> onepass;
    //vector<fpass> totalpass;
    fpass         temp;
    int pushnumber = 0;
    int k =0 ;
    int debug =0;
    vector<int> index;
    int beginCol;
    int u=0;
    vector<int> allRow;
    fpass  *passptr;
    elecol currentnode;
    char chainName;
    //cout<<"this the search program "<<endl;
/*
  if(a.size()<b.size())
  {
      ;
  }
  else
  //so the paramter a is whole elements in structure
  cout<<"the b size is "<<b.size()<<endl;
  for(i=0;i<b.size();i++)
     b[i].print();
  totalele[41].print();
  totalele[92].print();
  totalele[42].print();
  totalele[47].print();
*/
  {
    for(i=0;i<a.size();i++)
    {
        k = 0;
        if(a[i].geteleTyp() == b[k].geteleTyp()||b[k].geteleTyp()=='X')
        {

           currentnode = a[i];
           chainName = currentnode.getchain();
           passptr = new fpass(currentnode);
           passptr->addNodeToparent(currentnode);
           onepass.push(*passptr);
           delete passptr;
           while(!onepass.empty())
           {
              beginCol = onepass.front().getcurrentNode().getcol();
              debug++;
              k = onepass.front().getparent().size();
              for(j=beginCol+1;j<a.size();j++)
              {
                 if(a[j].geteleTyp() == b[k].geteleTyp()||b[k].geteleTyp()=='X')
                 {
                     currentnode = a[j];
                     //cout<<"the seconde chain name is "<<currentnode.getchain()<<endl;
                     if(chainjudge == true && chainName == currentnode.getchain())
                     {
                         passptr = new fpass(currentnode);
                         passptr->inheritPass(onepass.front());
                         passptr-> addNodeToparent(currentnode);
                         int y=0;
                         allRow.clear();
                         for(y=0;y<passptr->getparent().size();y++)
                         {
                            allRow.push_back(passptr->getparent()[y].getcol());
                            pushnumber++;
                         }
                         if(compareCol(allRow ,intM,row1,quM,row2,k) == true)
                         {
                            onepass.push(*passptr);
                         }
                         delete passptr;
                     }
                     if(chainjudge != true)
                     {
                         passptr = new fpass(currentnode);
                         passptr->inheritPass(onepass.front());
                         passptr-> addNodeToparent(currentnode);
                         int y=0;
                         allRow.clear();
                         for(y=0;y<passptr->getparent().size();y++)
                         {
                            allRow.push_back(passptr->getparent()[y].getcol());
                            pushnumber++;
                         }
                         if(compareCol(allRow ,intM,row1,quM,row2,k) == true)
                         {
                            onepass.push(*passptr);
                         }
                         delete passptr;
                     }
                 }
              }
              temp = onepass.front();
              if(temp.getparent().size()==row2)
                 totalpass.push_back(temp);
              onepass.pop();
           }
        }
    }
  }
    //cout<<"this out of search program"<<endl;
/*
    for(u=0;u<totalpass.size();u++)
    {
       cout<<"the no "<<u<<" times "<<endl;
       totalpass[u].print();
    }
*/
}   
void searchControl::printOuptfile(vector<fpass> &totalpass,vector<matrixElment> &IntMnumEl,char *pid,char *path)
{
   FILE *fileptr;
   vector<int> index;
   // The original 50-byte path1 overflows on hits-dir paths >49 chars (any
   // realistic NFS output path). The overflow then corrupts pid1/pid2
   // adjacent on the stack, and the subsequent strcat produces a garbled
   // filename that fopen silently fails on (or writes to the wrong place),
   // yielding "0 hits" sweep results even though the BFS found matches.
   // Sized for a typical NFS path + filename + slack.
   char path1[1024];
   char pid1[64];
   char pid2[200];
   if (strlen(path) >= sizeof(path1) - sizeof(pid1)) {
       cerr << "printOuptfile: hits-dir path too long ("
            << strlen(path) << " bytes), refusing to write: " << path << endl;
       return;
   }
   strncpy(path1, path, sizeof(path1) - 1);
   path1[sizeof(path1) - 1] = '\0';
   cout<<"this is printOuptfile function "<<endl;
   cout<<"the pid is "<<pid<<endl;
   cout<<"the path1 is "<<path1<<endl;
   int t = 0;
   //here I modify the file name, get the name end up with .ssd
   while(pid[t] != '\0')
     t++;
   pid[t-4] = '\0';
  // pid[t]='\0';
   if(int(totalpass.size())>0)
   {
      strcpy(pid1,"pdb");
      strcat(pid,".txt");
      strcat(pid1,pid);
      strcpy(pid2,pid1);
      t = 0;
      while(pid2[t] !='\0')
        t++;
      pid2[t-4]='\0';
      cout<<"after cat pid is "<<pid<<endl; 
      strcat(path1,pid1);
      fileptr = fopen(path1,"w");
      if(fileptr == NULL)
      {
         cout<<"the "<<pid1<<" can't open"<<endl;
         exit(0);
      }
      fprintf(fileptr,"%s\n",pid2);
      fprintf(fileptr,"sub-matrix:\n");
      int u;
      int i;
      int t;
      for(u=0;u<totalpass.size();u++)
      {
         fprintf(fileptr,"MOTIF:\n");
         for(i=0;i<totalpass[u].getparent().size();i++)
         {
            index.push_back(totalpass[u].getparent()[i].getcol());
         } 
         for(t=0;t<index.size();t++)
         {
             fprintf(fileptr,"segment-Type: %c Position:%4d Range:%5s--%5s %c Length:%4d\n",
                       IntMnumEl[index[t]].getEltype(),index[t]+1,(IntMnumEl[index[t]].getBeginNum()).c_str(),
                       (IntMnumEl[index[t]].getEndNum()).c_str(),IntMnumEl[index[t]].getChainName(),IntMnumEl[index[t]].getlenght());
         }
         index.clear();
      }
      fprintf(fileptr,"END\n");
      fclose(fileptr);
    }
}   
/*fprintf(fmatrix,"segment-Type: %c Position:%4d Range:%5d--%5d %c\n", segment[matc
hArray[i][l+1]].segType,matchArray[i][l+1]+1,segment[matchArray[i][l+1]].beginResId,segment[matchArray[i][l+1]].e
ndResID,segment[matchArray[i][l+1]].chainId);*/

bool searchControl::compareCol(vector<int> allRow,char **intM,int row1,char **quM,int row2, int nocol)
{
   int i;
   int col = allRow[allRow.size()-1];
   vector<char> onecolum;
   bool result = true;
   /*
   cout<<"the nocol is "<<nocol<<endl;
   cout<<"the col   is "<<col<<endl;
   cout<<endl;
   cout<<"the allRow is "<<endl;
   for(i=0;i<allRow.size();i++)
   cout<<allRow[i];
   cout<<endl;
   */
   for(i=0;i<allRow.size();i++)
   {
      onecolum.push_back(intM[allRow[i]][col]);
   }
   /*
   cout<<"onrcolum is "<<endl;
   for(i=0;i<onecolum.size();i++)
   cout<<onecolum[i]<<endl;
   for(i=0;i<allRow.size();i++)
   {
      cout<<allRow[i]<<" ";
   }
   cout<<endl;
   */
   if(nocol> row2)
   {
     cout<<"the qury matrix nocol parameter wrong "<<endl;
     exit(0);
   }
   /*
   if(onecolum.size() != nocol)
   {
     cout<<"onecolum size is "<<onecolum.size()<<endl;
     cout<<"this size maybe wrong "<<endl;
     exit(0);
   }
   */
   for(i=0;i<nocol;i++)
   {
      if(notequal(onecolum[i], quM[i][nocol]))
      {
          result = false;
          break;
      }
   }
   //cout<<"the result is "<<result<<endl;
   return result;
} 

//the char b means the letter of the qury matrix
//the char a means the letter of the interaction matrix
bool searchControl::notequal(char a , char b)
{
   bool result = true;
   if(a == b)
   {
      result = false;
      return result;
   }
   else
   {
      if(b == 'X')
        result = false;
      else if(b == 'x')
      {
        if(a == '-')
        {
          result = true;
        }
        else
          result = false;
      }
      else if(b == 'T')
      {
         if(a == 'v' || a == 't')
         {
             result = false;
         }
         if(a !='v' && a != 't')
         {
             result = true;
         }
      }
      else if(b == 'C')
      {
         if(a == 'u'||a == 'c')
         {
             result = false;
         }
         if(a !='u' && a!='c')
         {
             result = true;
         }
      }
      else
         result = true;
      return result;

    }
}       
#endif     
