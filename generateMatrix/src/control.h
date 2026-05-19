#ifndef CONTROL_H
#define CONTROL_H 
#include <iterator>
#include <algorithm>
#include <queue>
#include <map>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#if defined (BSD) && !_POSIX_SOURCE
    #include <sys/dir.h>
#else
    #include <dirent.h>
#endif
#include "element.h"
#include "sheet.h"
#include <dirent.h>
#include "mpi.h"
class control
{
    public :
           control();
           void readIndrefile(char *fileName , vector<element> &eleInstructure);
           void prodvector(vector<element> &eleInstru);
           void cpmyvector(myvector &a , myvector &b);//a -> b
           void printmyvector(myvector &a);
           void printPoint(point &a);
           void proInteractionMatr(vector<element> &eleInstru,FILE *matrixptr,char *pdbid);
           bool h_bond_E(element a , element b , bool& paralorNot);
           void prinLogfile(vector<element> &eleInstru, char *filename);
           void formsheet(vector<element> &eleInstru);
           void search(vector<element> &eleInstru , int i, vector<sheet> &sheetset,int sheetid);
           bool checkInsameSheet(int i , int j,vector<sheet> &sheetset);
           void Elmentdelete(vector<element> &eleInstru);
           void printInterMatrix(vector<element> &eleInstru , char **a,FILE *bug,char *pdbid);
           void prepareIndex(char *filename);
           void readIndrefile1(char *fileName , vector<element> &eleInstructure);
           void setElementHbond(vector<element> &eleInstructure,vector<h_bond> &pair,Fresidue *fptr);
           //these two functions are used to as sort the h_bond element
           void quicksort(vector<h_bond> &a,int start,int end);
           int  partition(vector<h_bond> &a ,int start,int end);
           void swap(vector<h_bond> &a,int left,int right);
           void bubblesort(vector<h_bond> &a);
           void manager(int numberprocess,char *option,char *filelist);
           bool formajob(char *filelist,char *sendarray,int &size);
           void worker(char *indrneelpathname,char *outputpathfile);
           int getjobnum(char *filelist,int filenumber);
           void setoption(char *option);
           void fillsheetdefault(vector<element> &ele);
           void sortElement(vector<element> &helix,vector<element> &strand,vector<element> &wholele);
           void bublesort(vector<element> &wholele);
           bool pairInmstrand(vector<element> &linker, h_bond *ptr);
           vector<element> groupEle(vector<element> &eleIn,int i);
           void resetSheetid(vector<element> &eleInstru, vector<element> &onesheet, vector<sheet> &sheetset,int &sheetid);
           ~control();
      private:
            vector<sheet> sheetsetdefault;
            map<string, int> pseudoAtomindex;
            map<int,string> rpseudoAtomindex;
            map<string,int> pseudoAtommark;
            char option[30];
};


control::control()
{
   ;
}
 
control::~control()
{
     ;
}

void control::fillsheetdefault(vector<element> &ele)
{
   int i=1;
   int j=0;
   sheet *sheettemp;
   bool judge = false;
   sheetsetdefault.clear(); 
   while(true)
   {
       judge = false;
       sheettemp = new sheet();
       for(j=0;j<ele.size();j++)
       {
          if(ele[j].getsheetId()==i)
          {
             sheettemp->addelement(j);
             judge = true;
          }
       } 
       i++;
       if(judge == true)
          sheetsetdefault.push_back(*sheettemp);
       delete sheettemp;
       if(judge == false)
         break;
    }
}   


void control::setoption(char *option1)
{
    strcpy(option,option1);
}


int control::getjobnum(char *filelist,int filenumber)
{
   int i = 0;
   int totalstar = 0;
   int jobnum =0;
   int startnum=0;
   cout<<"the filenuber is "<<filenumber<<endl;
  // cout<<filelist;
   while(filelist[i]!='\0')
   {
      if(filelist[i] == '*')
      {
         startnum++;
         totalstar++;
      }
      if(startnum == filenumber)
      {
         startnum =0 ;
         jobnum++;
      }
      i++;
   }
   if(startnum%filenumber !=0)
      jobnum++;
   cout<<"the totalstar is "<<totalstar<<endl;
   return jobnum; 
}
       
bool control::formajob(char *filelist,char *sendarray1,int &size)
{
    static int i=0;
    int j = 0;
    size = 0;
    int starnum = 0;
    while(filelist[i] !='\0')
    {
      sendarray1[j] = filelist[i];
      size++;
      if(filelist[i] == '*')
      {
         starnum++;
      }
      if(starnum == FILENUMBER)
      {
          sendarray1[j+1] = '\0';
          i++;
          break;
      }
      i++;
      j++;
    }
    sendarray1[size] = '\0';
    if(size>0)
       return true;
    else
       return false;
}
      
void control::manager(int numberprocess,char *option,char *filelist)
{
     char nameset[1000000];
     char tokenstring[20];
     int  count = 0;
     char sendarray[10000];
     int jobnum;
     int size;
     int startcount;
     int i,j;
     int sendjobnum=0;
     int dumb;
     int sender;
     MPI_Status status;
     dirent *listptr;
     DIR   *dp;
     if(strcmp(option,"-f")==0||strcmp(option,"-fs")==0)
     {
        ifstream inputfile(filelist,ios::in);
        if(!inputfile)
        {
           cout<<"the file "<<filelist<<" can't open"<<endl;
           MPI_Finalize();
           exit(0);
        }
        else
        {
           while (inputfile.getline(tokenstring ,20,'\n'))
           {
              strcat(tokenstring ,"*");
              strcat(nameset , tokenstring);
              count++;
           }           
        } 
     }
     if(strcmp(option,"-d")==0||strcmp(option,"-ds")==0)
     {
          //this I should read the directory file
         dp = opendir(filelist);
         if(dp == NULL)
         {
           cout<<"the file "<<filelist<<" can't open"<<endl;
           MPI_Finalize();
           exit(0);
         }
         else
         {
             for(;listptr=readdir(dp);)
             {
                if(strcmp(listptr->d_name,".") != 0 && strcmp(listptr->d_name,"..") != 0) 
                {
                    strcpy(tokenstring,listptr->d_name);
                    strcat(tokenstring ,"*");
                    strcat(nameset , tokenstring);
                }
             }
         }
     }
     jobnum = getjobnum(nameset,FILENUMBER);
     cout<<"the jobnumber is "<<jobnum<<endl;
     for(i=1;i<numberprocess;i++)
     {
        formajob((char*)nameset,(char*)sendarray,size);
       // cout<<"job is "<< sendarray<<"size is "<<size<<endl;
        MPI_Send(sendarray,SENDSIZE,MPI_CHAR,i,size,MPI_COMM_WORLD);
     }
     int w;
     for(w=0;w<jobnum;w++)
     {
        MPI_Recv( &dumb,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        sender = status.MPI_SOURCE;
        if(formajob((char*)nameset,(char*)sendarray,size) == true)
        {
            MPI_Send(sendarray,SENDSIZE,MPI_CHAR,sender,size,MPI_COMM_WORLD); 
        }
        else
        {
            MPI_Send(MPI_BOTTOM,0,MPI_DOUBLE,sender,0,MPI_COMM_WORLD);
        }
     }
     cout<<"process 0 out of loop"<<endl;
}
void control::worker(char *indrneelpathname,char *outpathnamefile)
{
    char getarray[10000];
    int myrank;
    int finish = 1;
    char filename[200];
    char tokenname[30];
    FILE *a;
    FILE *indrefile;
    int j = -1;
    int i=0;
    char command[100];
    vector<element> eleInOneStru;
    bool fin = false;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    sprintf(command,"%s%d",outpathnamefile,myrank);
    char filename1[200];
    strcpy(filename,indrneelpathname);
    a = fopen(command,"w");
    if(a == NULL)
    {
       cout<<"the file "<<command<<"can't open"<<endl;
       MPI_Finalize();
       exit(0);
    } 
    cout<<"the work rank is "<<myrank<<endl;
    MPI_Recv(getarray,SENDSIZE,MPI_CHAR,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    j=-1;
    for(i=0;i<status.MPI_TAG;i++)
    {
       j++;
       fin = false;
       tokenname[j] = getarray[i];
       if(tokenname[j] == '*')
       {
           tokenname[j] = '\0';
           fin = true;
           j = -1;
       }
       if(fin == true)
       {
          strcpy(filename1,filename);
          strcat(filename1,tokenname);
          eleInOneStru.clear();
          readIndrefile1(filename1 , eleInOneStru);
          Elmentdelete(eleInOneStru);
          prodvector(eleInOneStru);
          proInteractionMatr(eleInOneStru,a,tokenname);
       }
    }       
    cout<<"the MPI_ANY_TAG is "<<status.MPI_TAG<<endl;
    //cout<<"process "<<myrank<<" get "<<getarray<<endl;
    while(status.MPI_TAG > 0)
    {
       MPI_Send(&finish,1,MPI_INT,0,1,MPI_COMM_WORLD);
       MPI_Recv(getarray,SENDSIZE,MPI_CHAR,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
       j=-1;
       for(i=0;i<status.MPI_TAG;i++)
       {
         j++;
         fin = false;
         tokenname[j] = getarray[i];
         if(tokenname[j] == '*') 
         {       
            tokenname[j] = '\0'; 
            fin = true; 
            j = -1; 
         }       
         if(fin == true)
         {       
            strcpy(filename1,filename);
            strcat(filename1,tokenname);
            eleInOneStru.clear();
            readIndrefile1(filename1 , eleInOneStru);
            Elmentdelete(eleInOneStru);
            prodvector(eleInOneStru);
            proInteractionMatr(eleInOneStru,a,tokenname);
         }       
       }        
    }
    fclose(a);
    cout<<"process "<<myrank<<" out of loop "<<endl;
}        
    
//this function output is the process Indreneel file and get eleInstructure
//this data structure.
void control::Elmentdelete(vector<element> &eleInstru)
{
   int i;
   vector<element>::iterator p;
   p = eleInstru.begin();
   for(i=0;i<eleInstru.size();)
   {
       if(p->getresidueGroup().size()<2 && p->getElmentType() == "E")
       {
            eleInstru.erase(p);
       }
       else if(p->getresidueGroup().size()<5 && p->getElmentType() == "H")
       {
            //cout<< p->getbeginId()<<endl;
           
            eleInstru.erase(p);
       }
       else if(p->getresidueGroup().size()<=2 && p->getElmentType() == "L")
       {
            eleInstru.erase(p);
       }
       else
       {
           p = p + 1;
           i++;
       }
   }
   for(i=0;i<eleInstru.size();i++)
   {
       if(eleInstru[i].getresidueGroup().size()<2)
       {
           cout<<"in Elmentdelete function still contain one size element"<<endl;
           eleInstru[i].print();
           exit(0);
       }
   }
}
void control::prepareIndex(char *fileName)
{
   FILE *fevt;
   string wholeLine;
   char line[1000];
   string temp;
   string atomFlag;
   string hetatmFlag;
   int atomNum=1;
   string eleIndex;
   string oeleIndex;
   string chainName;
   string residueName;
   
   fevt = fopen(fileName,"r");
   if(fevt == NULL)
   {
       cout<<"the file "<<fileName<<" can't open"<<endl;
       MPI_Finalize();
       exit(0);
   }
   cout<<"prepare index"<<endl;
   while(fgets(line,1000,fevt)!=NULL)
   {
       //cout<<line;
       wholeLine = line;
    
       atomFlag=temp.assign(wholeLine,0,4);
       hetatmFlag=temp.assign(wholeLine,0,6);
    
       if(atomFlag == "ATOM" || hetatmFlag == "HETATM")
       {  
       	  temp.assign(wholeLine,12,4);
          if(temp == " CA ")
          {
             residueName = "CA";
             
             chainName = temp.assign(wholeLine,21,1);
             //cout<<"the chainName is "<<chainName<<endl;
             eleIndex=temp.assign(wholeLine,22,5);
             oeleIndex=eleIndex;
            
             eleIndex=eleIndex+chainName;
             //cout<<"the elementNum is "<<elementNum<<endl;
             if((pseudoAtomindex[eleIndex]) && (hetatmFlag == "HETATM")){
             		
             }else{	
            	 atomNum++;
             	 rpseudoAtomindex[atomNum]=oeleIndex;
               pseudoAtomindex[eleIndex]=atomNum;	
               if(atomFlag=="ATOM"){
               	 eleIndex=eleIndex+atomFlag;
               }
               if(hetatmFlag=="HETATM"){	
               	 eleIndex=eleIndex+hetatmFlag;
               }
               pseudoAtommark[eleIndex]=1;
             }
            //cout << "pseudoAtomindex " << pseudoAtomindex[eleIndex] << " " << eleIndex <<endl;
             
          }
       }     
       if(wholeLine.find("ENDMDL")==0 )
       { 
           break;
       }

   }
   fclose(fevt);
}
void control::readIndrefile1(char *fileName , vector<element> &eleInstructure)
{
   FILE *fevt;
   Fresidue *fptr;
   vector<element> hlix;
   vector<element> strand;
   vector<element> linker;
   vector<element> mstrand;
   element  *eptr;
   string wholeLine;
   char line[1000];
   string elementType;
   string chainName;
   string chainName1;
   string temp;
   int beginNum;
   string beginIndex;
   int endNum;
   string endIndex;
   int length;
   int elementNum;
   string pairIndex;
   string atomIndex;

   int sheetid;
   bool finish = false;
   bool judge = false;
   int sum;
   queue<h_bond> pair;
   vector<h_bond> pair1;
   char deleteD;
   int indexM;
   fevt = fopen(fileName,"r");
   if(fevt == NULL)
   {
       cout<<"the file "<<fileName<<" can't open"<<endl;
       MPI_Finalize();
       exit(0);
   }
   cout<<"this is the new file reader "<<endl;
   elementNum= -1;
   while(fgets(line,1000,fevt)!=NULL)
   {
       //cout<<line;
       wholeLine = line;
       if(wholeLine.find("HELIX")==0)
       {
           //temp.assign(wholeLine,0,5);
           elementType = "H";
           temp.assign(wholeLine,19,1);
           chainName = temp;
           eptr = new element(chainName,elementType);
           //cout<<"the temp name is "<<chainName; 
           beginIndex=temp.assign(wholeLine,21,5);
           beginIndex=beginIndex+chainName;
           //cout<<"the temp after 21 is "<<temp<<endl;
           beginNum = pseudoAtomindex[beginIndex];
           //cout<<"beginNum "<<beginNum;
           //<<endl;
           endIndex=temp.assign(wholeLine,33,5);
           endIndex=endIndex+chainName;
           endNum = pseudoAtomindex[endIndex];
           //cout<<"the endNUm  is "<<endNum<<endl;
           eptr->setStarEndnum(beginNum,endNum);
           temp.assign(wholeLine,71,5);
           //cout<<"temp is "<<temp<<endl;
           length = atoi(temp.c_str());
           // cout<<"the length is "<<length<<endl;
           eptr->setlength(length);
           hlix.push_back(*eptr);
           delete eptr;
       }
       if(wholeLine.find("SHEET")==0)
       {
           elementType = "E";
           temp.assign(wholeLine,21,1);
           /*
           cout<<"the temp is "<<temp<<endl; 
           if(temp == "m")
           {
              cout<<"the m happen "<<endl;
              exit(0);
           }
           */ 
           chainName = temp;
           eptr = new element(chainName,elementType);
           if(strcmp(option,"-f")==0||strcmp(option,"-ds")==0||strcmp(option,"-os")==0 || strcmp(option,"-d")==0
              || strcmp(option,"-o")==0)
           {
              temp.assign(wholeLine,11,3);
              //cout<<"the temp is "<<atoi(temp.c_str())<<endl;
              eptr->setsheetId(atoi(temp.c_str()));
              //cout<<"the id is "<<eptr->getsheetId()<<endl;
           }
             beginIndex=temp.assign(wholeLine,22,5);
           beginIndex=beginIndex+chainName;
           //cout<<"the temp after 21 is "<<temp<<endl;
           beginNum = pseudoAtomindex[beginIndex];
           //cout<<"beginNum "<<beginNum;
           //<<endl;
           endIndex=temp.assign(wholeLine,33,5);
           endIndex=endIndex+chainName;
           endNum = pseudoAtomindex[endIndex];
           eptr->setStarEndnum(beginNum,endNum);
           //this step I set up the m option
           temp.assign(wholeLine,70,1);
           deleteD = temp[0];
           eptr->setdeleteStr(deleteD);
           temp.assign(wholeLine,71,5);
           length = atoi(temp.c_str());
           // cout<<"the length is "<<length<<endl;
           eptr->setlength(length);
           strand.push_back(*eptr);
           if(deleteD == 'm')
              mstrand.push_back(*eptr);
           delete eptr;
        }
        if(wholeLine.find("LINKER") == 0)
        {
           //cout<<"this is linker section"<<endl;
           elementType = "L";
           temp.assign(wholeLine,19,1);
           chainName = temp;
           chainName1 = temp.assign(wholeLine,31,1);
           /*
           cout<<"the chainName is "<<chainName<<endl;
           cout<<"the chainName1 is "<<chainName1<<endl;
           if(chainName == "E")
           exit(0); 
           */
           if(chainName1 != chainName)
             elementNum--;
           else
           {
             eptr = new element(chainName,elementType);
             //cout<<"the temp name is "<<chainName<<endl;
             beginIndex=temp.assign(wholeLine,21,5);
             beginIndex=beginIndex+chainName;
             //cout<<"the temp after 21 is "<<temp<<endl;
             beginNum = pseudoAtomindex[beginIndex];
             //cout<<"beginNum "<<beginNum;
             //<<endl;
             endIndex=temp.assign(wholeLine,33,5);
             endIndex=endIndex+chainName;
             endNum = pseudoAtomindex[endIndex];
             eptr->setStarEndnum(beginNum,endNum);
             temp.assign(wholeLine,71,5);
             //cout<<"temp is "<<temp<<endl;
             length = atoi(temp.c_str());
             //cout<<"the length is "<<length<<endl;
             eptr->setlength(length);
             linker.push_back(*eptr);
             hlix.push_back(*eptr);
             delete eptr;
           }
        }
        if(finish==false && hlix.size() + strand.size() ==elementNum+1 && elementNum != -1)
        {
            finish = true;
            int h,s,e;
            h=0;s=0;
            sortElement(hlix,strand,eleInstructure);
        }
       //I put the linker element into the helix vector element
       if(wholeLine.find("HELIX")==0 || wholeLine.find("SHEET")==0 || wholeLine.find("LINKER") == 0)
          elementNum++;
       h_bond *hptr;
       int fit=0;
       double xcoord,ycoord,zcoord;
       Fresidue *fptr;
       h_bond *hptr1;
       string residueName;
       if(wholeLine.find("PAIRS")==0)
       {
           int countpair=0;
           int cindex;
           for(cindex=7;cindex<=12;cindex++)
           {
                if(wholeLine[cindex] != ' ')
                   countpair++;
           }
           //cout<<"the countpair is "<<countpair<<endl;
           if(countpair <= 3)
           {
               hptr = new h_bond();
               hptr1 = new h_bond();
               // cout<<"the wholeLine is "<<wholeLine<<endl;
               temp.assign(wholeLine,19,1);
               chainName = temp;
               // cout<<"this pair section"<<endl;
               //cout<<" the chainName is "<<chainName<<endl;
               pairIndex=temp.assign(wholeLine,21,5);
               pairIndex=pairIndex+chainName;
               elementNum = pseudoAtomindex[pairIndex];
               //cout<<"the elementNum is "<<elementNum<<endl;
               hptr->setMaster(chainName,elementNum);
               hptr1->setSuborder(chainName,elementNum);
               temp.assign(wholeLine,31,1);
               chainName = temp;
               // cout<<"the chainName is "<<chainName<<endl;
               pairIndex=temp.assign(wholeLine,33,5);
               pairIndex=pairIndex+chainName;
               elementNum =pseudoAtomindex[pairIndex];
               hptr->setSuborder(chainName,elementNum);
               hptr1->setMaster(chainName,elementNum);
               // cout<<"the second elementNum is "<<elementNum<<endl;
               // cout<<"the elementNum is "<<elementNum<<endl;
               // cout<<"the chainName is "<<chainName<<endl;

               if(pairInmstrand(mstrand,hptr)==false && pairInmstrand(mstrand,hptr1) == false)
               {
                  pair1.push_back(*hptr);
                  pair1.push_back(*hptr1);
               }
               delete hptr;
               delete hptr1;
            }
            if(countpair == 4)
            {
               hptr = new h_bond();
               hptr1 = new h_bond();
               //cout<<"the wholeLine is "<<wholeLine<<endl;
               temp.assign(wholeLine,20,1);
               chainName = temp;
               //cout<<"this pair section"<<endl;
               //cout<<" the chainName is "<<chainName<<endl;
               pairIndex=temp.assign(wholeLine,22,5);
               pairIndex=pairIndex+chainName;
               elementNum = pseudoAtomindex[pairIndex];
               hptr->setMaster(chainName,elementNum);
               hptr1->setSuborder(chainName,elementNum);
               temp.assign(wholeLine,32,1);
               chainName = temp;
               // cout<<"the chainName is "<<chainName<<endl;
               pairIndex=temp.assign(wholeLine,34,5);
               pairIndex=pairIndex+chainName;
               elementNum = pseudoAtomindex[pairIndex];
               hptr->setSuborder(chainName,elementNum);
               hptr1->setMaster(chainName,elementNum);
               // cout<<"the second elementNum is "<<elementNum<<endl;
               // cout<<"the elementNum is "<<elementNum<<endl;
               // cout<<"the chainName is "<<chainName<<endl;
               if(pairInmstrand(mstrand,hptr)==false && pairInmstrand(mstrand,hptr1) == false)
               {
                  pair1.push_back(*hptr);
                  pair1.push_back(*hptr1);
               }
               //pair1.push_back(*hptr);
               //pair1.push_back(*hptr1);
               delete hptr;
               delete hptr1;
            }
            if(countpair > 4)
            {
               cout<<"in indranel1 file reader is wrong "<<endl;
               exit(0);
            }
       }
       if(wholeLine.find("ATOM")==0 || wholeLine.find("HETATM")==0)
       {
          if(judge == false && (int)pair1.size()>0)
          {
              bubblesort(pair1);
              /*
              for(int w=0;w<pair1.size();w++)
                 pair1[w].print();
              exit(0);
              */
              judge = true;
          
              for(int y=0;y<pair1.size();y++)
              {
                  pair.push(pair1[y]);
              }
          }
          temp.assign(wholeLine,12,4);
          if(temp == " CA ")
          {
             residueName = "CA";
             chainName = temp.assign(wholeLine,21,1);
             //cout<<"the chainName is "<<chainName<<endl;
             atomIndex=temp.assign(wholeLine,22,5);
             atomIndex=atomIndex+chainName;
             elementNum = pseudoAtomindex[atomIndex];
             if(wholeLine.find("ATOM")==0){
             	atomIndex=atomIndex+"ATOM";
             }	
             if(wholeLine.find("HETATM")==0){
             	atomIndex=atomIndex+"HETATM";
             }	
             if(pseudoAtommark[atomIndex]){
             	
             }else{
             	 continue;
             }
             //cout<<"the elementNum is "<<elementNum<<endl;
             temp.assign(wholeLine,30,8);
             xcoord = atof(temp.c_str());
             //cout<<"the xcoord is "<<xcoord<<endl;
             temp.assign(wholeLine,38,8);
             ycoord = atof(temp.c_str());
             //cout<<"the ycoord is "<<ycoord<<endl;
             temp.assign(wholeLine,46,8);
             zcoord = atof(temp.c_str());
             //cout<<"the zcoord is "<<zcoord<<endl;
             fptr = new Fresidue(chainName , elementNum,residueName,xcoord ,ycoord,zcoord," ");
             setElementHbond(eleInstructure,pair1,fptr);
             //cout<<"after setElementHbond(eleInstructure,pair,fptr) function"<<endl;
             delete fptr;
          }
       }     
       if(wholeLine.find("ENDMDL")==0 )
       { 
           break;
       }

   }
   fclose(fevt);
}

//this funciton I will junk if ptr is belong to m strand
bool control::pairInmstrand(vector<element> &mstrand, h_bond *ptr)
{
  int i,j;
  int renum;
  char chid;
  bool find = false;
  for(i=0;i<mstrand.size();i++)
  {
     if(mstrand[i].withelement(ptr->getMasterNum(),ptr->getMasterChain()[0]) == true)
     {
        find = true;
        break;
     }
  }
  return find;
}

void control::sortElement(vector<element> &hlix,vector<element> &strand,vector<element> &eleInstructure)
{
   int h,s,e,i;
   vector<char> chainName;
   vector<element> hlix1;
   vector<element> wholelement;
   char id;
   bool first = true;
   bool find = false;

   hlix1 = hlix;
   for(h=0;h<strand.size();h++)
   {
      hlix1.push_back(strand[h]);
   }
   //this code I deal with the chainId about the helix
   for(h=0;h<hlix1.size();h++)
   {
      id = hlix1[h].getChainId()[0];
      if(first== true) 
      {
         chainName.push_back(id);
         first = false;
      }
      find = false;
      for(s=0;s<chainName.size();s++)
      {
         if(chainName[s] == id)
         {
            find = true;
            break;
         }
      }
      if(find == false)
         chainName.push_back(id);
   } 
   //Now chainName already has all the record of the different chain names
   if(chainName.size() == 0)
   {
      cout<<"in sortElement function, the chainName status is wrong"<<endl;
      exit(0);
   }
   /*
   cout<<"chainName size is "<<chainName.size()<<endl;
   
   for(h=0;h<chainName.size();h++)
      cout<<chainName[h]<<endl;
   */
   for(h=0;h<chainName.size();h++)
   {
     id = chainName[h]; 
     wholelement.clear();
     for(s=0;s<hlix1.size();s++)
     {
       if(hlix1[s].getChainId()[0] == id)
       {
          wholelement.push_back(hlix1[s]); 
       }
     }
     bublesort(wholelement);
     for(i=0;i<wholelement.size();i++)
        eleInstructure.push_back(wholelement[i]);
   }
     
     
     
}
void control::bublesort(vector<element> &a)
{
   element hold;
   int i;
   for(int pass=0;pass<a.size()-1;pass++)
      for(i=0;i<a.size()-1;i++)
      {
         if(a[i].getbeginNum()>a[i+1].getbeginNum())
         {
             hold = a[i];
             a[i] = a[i+1];
             a[i+1] =hold;
         }
       }
}

void control::bubblesort( vector<h_bond>&a)
{
   h_bond hold;
   int i;
   for(int pass=0;pass<a.size()-1;pass++)
      for(i=0;i<a.size()-1;i++)
      {
         if(a[i]>a[i+1])
         {
             hold = a[i];
             a[i] = a[i+1];
             a[i+1] =hold;
         }
       }
}

void control::swap(vector<h_bond> &a,int left,int right)
{
   h_bond temp;
   temp = a[left];
   a[left] = a[right];
   a[right] = temp;
}
 
//this is function I put the residue into the element and at the same time, I add the h_bond 
//relation into the residue
void control::setElementHbond(vector<element> &eleInstructure,vector<h_bond> &pair,Fresidue *fptr)
{
   int i;
   int master;
   string maschain;
   int suborder;
   string subchain;
   /*
   while(!pair.empty())
   {
      if(pair.front().getMasterChain() == fptr->getchaiName() &&
         pair.front().getMasterNum() == fptr->getRsidNumber())
       {
          fptr->addHbondResidue(pair.front());
          pair.pop();
       }
      else     
         break;
    }
    */
   for(i=0;i<pair.size();i++)
   {
      if(pair[i].getMasterChain() == fptr->getchaiName() && pair[i].getMasterNum() == fptr->getRsidNumber())
      {
         fptr->addHbondResidue(pair[i]);
      }
   }
   for(i=0;i<eleInstructure.size();i++)
   {
      if(fptr->getRsidNumber()>=eleInstructure[i].getbeginNum() && fptr->getRsidNumber()<=eleInstructure[i].getendNum()&&fptr->getchaiName()==eleInstructure[i].getChainId())
      {
         if(eleInstructure[i].getopresidue() != fptr->getRsidNumber())
         eleInstructure[i].addResidue(*fptr);
      }
   }
} 

void control::readIndrefile(char *fileName , vector<element> &eleInstructure)
{  
   FILE *fevt;
   Fresidue *ptr;
   element  *eptr=NULL;
   bool con = true;
   char oneLine[200];
   string wholestring;
   bool   samelment = false;
   string initialElment = "z";
   int atomNum;
   string cname;
   int Rnum;
   string Rname;
   bool smallv = false;
   bool bigO = false;
   double x,y,z, disiToi3 , Torsion ;
   string eleNam = " " ;
   h_bond *hptr; 
   string temp;
   int i ;
   int count = 0; 
   string prestring;
   string H_cname;
   int    H_Rnum;
   fevt = fopen(fileName ,"r");
   int tabposition;
   if(fevt == NULL)
   {
       cout<<"in control class read indrel file can't open "<<endl;
       exit(0);
   }
   else
   { 
    while(con == true)
    {
     if(fgets(oneLine ,200 ,fevt) != NULL)
     {
       int tabsit = 0;
       wholestring=oneLine;
       cname.assign(wholestring , 7 ,1); 
       if(cname == " ")
       {
          cname = "A";
       }
       i=0;
     
       while(oneLine[i] !='\0')
       {
           if(oneLine[i] == '\t')
           {
               tabposition = i;
               break;
           }
           i++;
       }
      
       temp.assign(wholestring,8 , 5);
       Rnum = atoi(temp.c_str());
       Rname.assign(wholestring , 15 ,3);
       temp.assign(wholestring ,19 ,8);
       x = atof(temp.c_str());
       temp.assign(wholestring,28 ,8);
       y = atof(temp.c_str());
       temp.assign(wholestring ,37 ,8);
       z = atof(temp.c_str());
       temp.assign(wholestring ,tabposition+17,1);
       //cout<<"wholeline size is "<<wholestring.size();
       //cout<<wholestring<<endl;
          if(temp == "O" || temp == "v")
          {
             prestring = eleNam;
          }
          eleNam = temp;
          if(eleNam == " ")
          {
              eleNam = "s"; 
          }
          
       //   if(Rnum == 35 &&eleNam == "s")
       //   exit(0);
          if(wholestring.size() == tabposition+17+2)
          {
             ptr = new Fresidue(cname ,Rnum, Rname,x,y,z,eleNam); 
          }
          if(wholestring.size() == tabposition+26)
          {
              H_cname.assign(wholestring,tabposition+19,1);
              if(H_cname == " ")
                 H_cname = "A";
              temp.assign(wholestring , tabposition+20,5);
              H_Rnum = atoi(temp.c_str());
            
              hptr = new h_bond(H_cname , H_Rnum);
              ptr = new Fresidue(cname ,Rnum, Rname,x,y,z,eleNam);
              ptr->addHbondResidue(*hptr);
              delete hptr;
          }
          if(wholestring.size() == tabposition+33)
          {
             H_cname.assign(wholestring,tabposition+19,1);
             if(H_cname == " ")
                 H_cname = "A";
             temp.assign(wholestring ,tabposition+20 ,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr = new Fresidue(cname ,Rnum, Rname,x,y,z,eleNam);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             H_cname.assign(wholestring , tabposition+26 ,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+27,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
          }
          if(wholestring.size() == tabposition + 40)
          {
             H_cname.assign(wholestring,tabposition+19,1);
             if(H_cname == " ")
                 H_cname = "A";
             temp.assign(wholestring ,tabposition+20 ,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr = new Fresidue(cname ,Rnum, Rname,x,y,z,eleNam);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             H_cname.assign(wholestring , tabposition+26 ,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+27,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             H_cname.assign(wholestring,tabposition+33,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+34,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
          }
          if(wholestring.size() > 109 )
          {
             H_cname.assign(wholestring,tabposition+19,1);
             if(H_cname == " ")
                 H_cname = "A";
             temp.assign(wholestring ,tabposition+20 ,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr = new Fresidue(cname ,Rnum, Rname,x,y,z,eleNam);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             H_cname.assign(wholestring , tabposition+26 ,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+27,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             H_cname.assign(wholestring,tabposition+33,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+34,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             
             H_cname.assign(wholestring,tabposition+40,1);
             if(H_cname == " ")
             {
                H_cname = "A";
             }
             temp.assign(wholestring,tabposition+41,5);
             H_Rnum = atoi(temp.c_str());
             hptr = new h_bond(H_cname , H_Rnum);
             ptr->addHbondResidue(*hptr);
             delete hptr;
             cout<<"the H_cname is"<<H_cname<<endl;
             cout<<"the H_Rnum is "<<H_Rnum<<endl;
             
          }
          if(wholestring.size()>116)
          {
             
             cout<<"the tabposition is "<<tabposition<<endl; 
             cout<<"the wholestring size is "<<wholestring.size()<<endl;
             cout<<"in readindreneel file wrong wholestring size wrong "<<endl;
             exit(0);
          }
          bool smallo = false;
          if(initialElment != eleNam )
          {
             if(eleNam == "s")
             {
                 if(eptr != NULL)
                 {
                     eleInstructure.push_back(*eptr);
                     delete eptr;
                     eptr = NULL;
                     initialElment = "z";
                 }
                 prestring = "z";
                 smallv = false;
             }
             if(eleNam == "o")
             {
                  if(eptr != NULL)
                  {
                     ptr->setEletyp("H");
                     eptr->setEltype("H");
                     eptr->addResidue(*ptr);
                     eleInstructure.push_back(*eptr);
                     delete eptr;
                     eptr = new element(cname , "H");
                     eptr->addResidue(*ptr);
                     delete ptr;
                     initialElment = "H";
                     smallo = true;
                  }
                  else
                  { 
                      cout<<"in small o section wrong "<<endl;
                      exit(0);
                  }
                 prestring = "z";
                 smallv = false;
             }
            if(eleNam == "n")
            {
                ptr->setEletyp("H");
             //   eptr->setEltype("H");
                eleInstructure[eleInstructure.size()-1].addResidue(*ptr);
                if(eptr == NULL)
                {
                    cout<<"in n section something is wrong "<<endl;
                    exit(0);
                }
                else
                {  
                    ptr->setEletyp("E");
                    eptr->addResidue(*ptr); 
                    
                }
                delete ptr;
             }
            if(eleNam == "v"&& prestring != "O" )
            {
                element *eptr1;
                ptr->setEletyp("H");
                eptr->addResidue(*ptr);
                eleInstructure.push_back(*eptr);
                delete eptr; 
                eptr1 = new element(cname,"E");
                eptr1->addResidue(*ptr);
                eleInstructure.push_back(*eptr1);
                delete eptr1;
                eptr = new element(cname,"H");
                eptr->addResidue(*ptr);
//              prestring = "O";
            }
            if(eleNam == "v" && prestring == "O")
            {
                cout<<" v and prestring == O happen "<<endl;
                ptr->setEletyp("H");
                eleInstructure[eleInstructure.size()-1].addResidue(*ptr); 
                ptr->setEletyp("E");
                eptr->addResidue(*ptr);
                eleInstructure.push_back(*eptr);
                delete eptr;
                eptr = new element(cname,"H");
                ptr->setEletyp("H");
                eptr->addResidue(*ptr);
                delete ptr;
            } 
            if(eleNam == "O")
            {
                if(eptr != NULL && prestring != "O" && prestring != "n"&&prestring != "o"&&
                   prestring != "v")
                {
                   ptr->setEletyp(prestring);
                   eptr->setEltype(prestring); 
                   eptr->addResidue(*ptr);
                   eleInstructure.push_back(*eptr);
                   delete eptr;
                   eptr = new element(cname , prestring);
                   if(prestring == "H")
                   {
                        ptr->setEletyp("E");
                        eptr->setEltype("E");
                        initialElment = "E";
                   }
                   else if(prestring == "E")
                   {
                        ptr->setEletyp("H");
                        eptr->setEltype("H");
                        initialElment = "H";
                   }
                   else 
                   {
                        cout<<"prestring state is wrong "<<endl;
                        cout<<"the current state of prestring is "<<prestring<<endl;
                        exit(0);
                   }
                   eptr->addResidue(*ptr);
                   delete ptr;
                 }
                 if(eptr == NULL)
                 {
                     if(prestring == "V")
                     {
                        eptr = new element(cname , "H");
                        ptr->setEletyp("H");
                        eptr->addResidue(*ptr);
                        initialElment = "H";
                     }
                     else if(prestring == "s")
                     {
                         element *eptr1 = new element(cname , eleNam);
                         eptr1->addResidue(*ptr);
                         eleInstructure.push_back(*eptr1);
                         delete eptr1;
                         eptr = new element(cname , eleNam); 
                         eptr->addResidue(*ptr);
                         bigO = true;
                         initialElment = "z"; 
                         delete ptr;
                     }
                     else if(prestring == " " || prestring == "e")
                     {
                         element *eptr1 = new element(cname , eleNam);
                         eptr1->addResidue(*ptr);
                         eleInstructure.push_back(*eptr1);
                         delete eptr1;
                         eptr = new element(cname , eleNam); 
                         eptr->addResidue(*ptr);
                         bigO = true;
                         initialElment = "z"; 
                         delete ptr;
                     }
                     else 
                     {
                         
                         cout<<"in big O, eptr == NULL, prestring is wrong "<<endl;
                         cout<<"prestring is "<<prestring<<endl;
                         cout<<"in big O section wrong  "<<endl;
                         exit(0);
                      }
                 }
          
                 if(prestring == "O")
                 {
                   if(smallv == false) 
                   {
                    if(initialElment == "E")
                       ptr->setEletyp("H");
                    else if(initialElment == "H")
                       ptr->setEletyp("E");
                    else if(initialElment == "z")
                    {
                      if(bigO == false) 
                        initialElment = "H";
                      else
                      {
                        initialElment = "z";
                      }
                    }
                    else 
                    {
                       cout<<"initail wrong "<<endl;
                       cout<<"the initialELment is "<<initialElment<<endl;
                       exit(0);
                    }
                    
                   } 
                  
                   else
                   {
                          ;
                   }
                    eleInstructure[eleInstructure.size()-1].addResidue(*ptr); 
                    eptr->addResidue(*ptr);
                    eptr->setEltype(initialElment);
                    delete ptr;
                 }
   
                 if(prestring == "n")
                 {
                    ptr->setEletyp("E");
                    eptr->addResidue(*ptr);
                    eleInstructure.push_back(*eptr);
                    delete eptr;
                    eptr = new element(cname , "H");
                    ptr->setEletyp("E");
                    eptr->addResidue(*ptr);
                    initialElment = "H";
                    delete ptr;
                 }
                 if(prestring == "o" && eptr != NULL)
                 {
                    element *eptr2 = new element(cname, eleNam);
                    eptr2->addResidue(*ptr);
                    eptr2->setEltype("E");
                    eleInstructure.push_back(*eptr2);
                    delete eptr2;
                    eptr->addResidue(*ptr);
                    initialElment = "H";
                 }   
                 if(prestring == "v" && eptr != NULL)
                 {
                    ptr->setEletyp("E");
                    eleInstructure[eleInstructure.size()-1].addResidue(*ptr);
                    ptr->setEletyp("H");
                    eptr->addResidue(*ptr);
                    initialElment = "H";
                    smallv = true;
                  }
                    
                    
                    
            } 
            if(eleNam == "e")
            {
                if(eptr != NULL)
                {
                    ptr->setEletyp("H");
                    eptr->addResidue(*ptr);
                    eleInstructure.push_back(*eptr);
                    delete eptr;
                    eptr = NULL;
                    initialElment = "z";
                }
                else
                {
                    cout<<"in small e section wrong  "<<endl;
                    exit(0);
                }
                prestring = "z";
                smallv = false; 
             }
            if(eleNam == "V")
            {
                if(eptr != NULL)
                {
                    ptr->setEletyp("E");
                    eptr->addResidue(*ptr);
                    eleInstructure.push_back(*eptr);
                    initialElment = "z";
                    delete eptr;
                    eptr = NULL;
                }         
                else
                {
                    cout<<"in big V section   wrong  "<<endl;
                    exit(0);
                }
                prestring = "z";
                smallv = false;
             }
            if(eleNam == "H")
            {
/*
              if(Rnum == 41)
              {
                 cout<<"the bigO is "<<bigO<<endl;
                 exit(0);
              }
 */             
              if(bigO == true)
              {
                 eleInstructure[eleInstructure.size()-1].setEltype("E");
                 eptr->setEltype("H");
                 eptr->addResidue(*ptr);
                 delete ptr;
                 initialElment = "H";
                 bigO = false;
              }
              else
              {
                if(eptr != NULL)
                {
                    eleInstructure.push_back(*eptr);
                    if(smallo == false)
                    delete eptr;
                    if(smallo == false)
                    eptr = new element(cname ,eleNam);
                    eptr->addResidue(*ptr);
                    initialElment = "H";
                }
                else
                {
                   
                    if(Rnum == 41 && eptr ==NULL)
                    {
                      cout<<"the smallo is "<<smallo<<endl;
                    }
                    eptr = new element(cname ,eleNam);
                    eptr->addResidue(*ptr);
                    initialElment = "H";
                    delete ptr;
                }
                if(smallo == true)
                    smallo = false;
              }
              prestring = "z";
              smallv = false;
             }
            if(eleNam == "E")
            {
              if(bigO ==true)
              {
                 eleInstructure[eleInstructure.size()-1].setEltype("H");
                 eptr->setEltype("E");
                 eptr->addResidue(*ptr);
                 delete ptr;
                 initialElment = "E";
                 bigO = false;
              }
              else
              {
               if(eptr != NULL)
               {
                   eleInstructure.push_back(*eptr);
                   delete eptr;
                   eptr = new element(cname ,eleNam);
                   eptr->addResidue(*ptr);
                   initialElment = "E";
               }
               else
               { 
                    eptr = new element(cname ,eleNam);
                    eptr->addResidue(*ptr);
                    initialElment = "E";
               }
              }
              prestring = "z";
              smallv = false;
             }
             if(eleNam == " ")
             {
                if(eptr !=NULL)
                {
                     eleInstructure.push_back(*eptr);
                     delete eptr;
                     eptr = NULL;
                }
                initialElment = "z";
                prestring = "z";
                smallv = false;
             }
          }
         else
         {
               if(eleNam != "s"&&eleNam !=" ")
                  eptr->addResidue(*ptr);
         }
       }
     else
      {
         con = false;
         if(eptr != NULL)
         {
             eleInstructure.push_back(*eptr);
             delete eptr;
         }
         fclose(fevt);
         
      }
    }
     
  }
}


bool control::checkInsameSheet(int i , int j,vector<sheet> &sheetset)
{
    int k=0;
    int m=0;
    int size;
    bool insame1 = false;
    bool insame2 = false;
    for(k=0 ;k<sheetset.size();k++)
    {
        insame1 = false;
        insame2 = false; 
        size = sheetset[k].getSheet().size();
        for(m=0;m<size;m++)
        {
            if(i == sheetset[k].getSheet()[m])
               insame1 = true;
            if(j == sheetset[k].getSheet()[m])
               insame2 = true;
         }
        if(insame1 == true && insame2 == true)
           break;
    }
        if(insame1 == false && insame2 == false)
            return false;
        if(insame1 == false && insame2 == true)
            return false;
         if(insame1 == true && insame2 == false)
            return false;
        if(insame1 == true && insame2 == true)
            return true;
}    
//this function output is I use the vector to represent the elements.
void control::prodvector(vector<element> &eleInstru)
{
   int i , j;
   double *all_H_cordinates;
   point *projectarray;
   myvector parameter;
   int size;
   int reIndex =0;
   int k =0;
   int m;
   for(i=0 ; i<eleInstru.size();i++)
   {
         size = eleInstru[i].getResidueSize();
         reIndex = 0;
         all_H_cordinates = new double[3*size];
         projectarray     = new point[size];
         for(j=0 ;j<size; )
         {
              all_H_cordinates[reIndex] = eleInstru[i].getresidueGroup()[j].get_X();
               //  cout<<eleInstru[i].getresidueGroup()[j].get_X()<<endl;
              reIndex++;
              all_H_cordinates[reIndex] = eleInstru[i].getresidueGroup()[j].get_Y();
               //  cout<<eleInstru[i].getresidueGroup()[j].get_Y()<<endl;
              reIndex++;
              all_H_cordinates[reIndex] = eleInstru[i].getresidueGroup()[j].get_Z();
              //   cout<<eleInstru[i].getresidueGroup()[j].get_Z()<<endl;
              reIndex++;
              j++;
         }
         if(eleInstru[i].getElmentType() == "H" && size >3)
         {
             rotFit(all_H_cordinates ,size ,3, &parameter);
             eleInstru[i].setVector(parameter);
             cpmyvector(eleInstru[i].getVector() ,parameter); 
             //printmyvector(parameter);
             pointProjToVect(all_H_cordinates,size,parameter,projectarray);
         }
         if(eleInstru[i].getElmentType() == "E" && size > 3 || eleInstru[i].getElmentType() == "L" && size > 3)
         {
             k = (int) (size/3);
             parameter.initPoint.xcoord = (all_H_cordinates[(k-1)*3 + 0]+all_H_cordinates[k*3 +0])/2;
             parameter.initPoint.ycoord = (all_H_cordinates[(k-1)*3 + 1]+all_H_cordinates[k*3 +1])/2; 
             parameter.initPoint.zcoord = (all_H_cordinates[(k-1)*3 + 2]+all_H_cordinates[k*3 +2])/2;
             parameter.endPoint.xcoord = (all_H_cordinates[(size-k)*3 +0] +all_H_cordinates[(size-k-1)*3 +0])/2;
             parameter.endPoint.ycoord = (all_H_cordinates[(size-k)*3 +1] +all_H_cordinates[(size-k-1)*3 +1])/2;
             parameter.endPoint.zcoord = (all_H_cordinates[(size-k)*3 +2] +all_H_cordinates[(size-k-1)*3 +2])/2;
             eleInstru[i].setVector(parameter);
             pointProjToVect(all_H_cordinates,size,parameter,projectarray);
             parameter.initPoint.xcoord = projectarray[0].xcoord;
             parameter.initPoint.ycoord = projectarray[0].ycoord;
             parameter.initPoint.zcoord = projectarray[0].zcoord;
             parameter.endPoint.xcoord  = projectarray[size-1].xcoord;
             parameter.endPoint.ycoord  = projectarray[size-1].ycoord;
             parameter.endPoint.zcoord  = projectarray[size-1].zcoord;
             //printmyvector(parameter);
             eleInstru[i].setVector(parameter);
         }
         if(eleInstru[i].getElmentType() == "E" && size == 3 || eleInstru[i].getElmentType() == "L" && size == 3) 
         {
             parameter.initPoint.xcoord = (all_H_cordinates[0] +all_H_cordinates[3])/2;
             parameter.initPoint.ycoord = (all_H_cordinates[1] +all_H_cordinates[4])/2;
             parameter.initPoint.zcoord = (all_H_cordinates[2] +all_H_cordinates[5])/2;
             parameter.endPoint.xcoord = (all_H_cordinates[3] +all_H_cordinates[6])/2;
             parameter.endPoint.ycoord = (all_H_cordinates[4] +all_H_cordinates[7])/2;
             parameter.endPoint.zcoord = (all_H_cordinates[5] +all_H_cordinates[8])/2;
             pointProjToVect(all_H_cordinates,size,parameter,projectarray);
             parameter.initPoint.xcoord = projectarray[0].xcoord;
             parameter.initPoint.ycoord = projectarray[0].ycoord;
             parameter.initPoint.zcoord = projectarray[0].zcoord;
             parameter.endPoint.xcoord = projectarray[2].xcoord;
             parameter.endPoint.ycoord = projectarray[2].ycoord;
             parameter.endPoint.zcoord  = projectarray[2].zcoord;
             //printmyvector(parameter);
             eleInstru[i].setVector(parameter);
          } 
         if(eleInstru[i].getElmentType() == "E" && size == 2)
         {
             parameter.initPoint.xcoord = all_H_cordinates[0];
             parameter.initPoint.ycoord = all_H_cordinates[1];
             parameter.initPoint.zcoord = all_H_cordinates[2];
             parameter.endPoint.xcoord = all_H_cordinates[3];
             parameter.endPoint.ycoord = all_H_cordinates[4];
             parameter.endPoint.zcoord = all_H_cordinates[5];
             projectarray[0] = parameter.initPoint;
             projectarray[1] = parameter.endPoint;
         }
          
        for(m=0 ;m<size;m++)
        {
            eleInstru[i].addProPoint(projectarray[m]);
        }
        delete [] all_H_cordinates;
        delete [] projectarray;
   }
    
} 
      
//this function I will get the elment that belong to same sheet,sheet id is i
//parameter i is the sheetid, also In this function the return elements I will get rid of the strand with m mark
vector<element> control::groupEle(vector<element> &eleIn,int i)
{
   vector<element> oneSheet;
   vector<int> elementId;
   int j;

   elementId.clear();
   oneSheet.clear();
   //sheetsetdefault is a globe vairable,vector<sheet> which have element index in each of sheet structure
   
   elementId = sheetsetdefault[i].getSheet();
   for(j=0;j<elementId.size();j++)
   {
      if(eleIn[elementId[j]].getdeleteStr() != 'm')
      {
         oneSheet.push_back(eleIn[elementId[j]]);
      }
      else
        eleIn[elementId[j]].setsheetId(0);
      /*
      if(eleIn[elementId[j]].getdeleteStr() == 'm')
         eleIn[elementId[j]].print();
      */
   }
   return oneSheet;
}

//this function I will remap the sheetid into eleInstru data,if strand in elestru is m sheetid is 0, so is 
//sheet with one strand
void control::resetSheetid(vector<element> &eleInstru, vector<element> &onesheet, vector<sheet> &sheetset,int &sheetid)
{
   int i,j;
   int eleIndex;
   for(i=0;i<sheetset.size();i++)
   {
     if(sheetset[i].getSheet().size() == 1)
     {
        eleIndex = sheetset[i].getSheet()[0]; 
        onesheet[eleIndex].setsheetId(0);
        sheetid--;
     }
   }
   for(i=0;i<onesheet.size();i++)
   {
      for(j=0;j<eleInstru.size();j++)
      {
         if(onesheet[i] == eleInstru[j])
            eleInstru[j].setsheetId(onesheet[i].getsheetId());
      }
   }
   /*
   for(i=0;i<onesheet.size();i++)
     cout<<onesheet[i].getsheetId()<<" ";
   cout<<endl;
   exit(0);
   */
}

//this function output is a matrix which desxribe the relationship
//between elements.
void control::proInteractionMatr(vector<element> &eleInstru, FILE *matrixptr, char * pdbid)
{
    int i,j,k;
    char bugfile[20];
    FILE *bug;
    int size,size1,size2;
    char **intMatrix;
    double overlap;
    double distance;
    double anglevalue;
    myvector v1 , v2;
    point p1 ,p2 ,p3 ,p4;
    string trueBeginpos;
    string trueEndpos;
    //char trueBeginpos[100];
   	//char trueEndpos[100];
    int pseudoBeginpos;
    int pseudoEndpos;
    int thiseleLength;
    string begin1;
    string end1;
    string chainB1;
    string chainE1;
    vector<sheet> sheetset;
    int indexS;
    point *elecoord, *elecoord2;
    bool paralorNot = false;
    static FILE *matrixfile;
    vector<Fresidue> a1 , b1;
    vector<h_bond> H_bond1;
    vector<element> onesheetEle;
    int sheetid = 0;
    int idindex,w;

    size = eleInstru.size();
    cout<<"the size is "<<size<<endl;
    intMatrix = new char *[size];
    for(i=0 ; i<size;i++)
    {
       intMatrix[i] = new char [size];
    }
    cout<<"the eleInstru size is "<<eleInstru.size()<<endl;
    //exit(0);
    //this step I will generate the sheet informatin by myself
    if(strcmp(option ,"-o")==0||strcmp(option ,"-f")==0||strcmp(option ,"-d")==0)
    //if(strcmp(option ,"-o")==0||strcmp(option ,"-f")==0||strcmp(option ,"-d")==0)	originally, -d write as -ds.
    {
       //this is step I initialize the sheetsetdefault data using indraneel files
       fillsheetdefault(eleInstru);
       //cout<<"the sheetsetdefault size is "<<sheetsetdefault.size()<<endl;
       //vector<sheet> sheetsetdefault is the vector sheet data structure
       for(indexS=0;indexS<sheetsetdefault.size();indexS++)
       {
          onesheetEle.clear();

          //this function I return a group of element with the same sheet Id
          //sheetsetdefault is a globe vairable,vector<sheet> which have element index in each of sheet structure
          //groupEle function access this varible

          onesheetEle = groupEle(eleInstru,indexS);
          /*
          if(indexS == 1)
          {
             cout<<"onesheetEle size is "<<onesheetEle.size()<<endl;
             for( w=0;w<onesheetEle.size();w++)
                onesheetEle[w].print();
             exit(0);
          }
          */
          //onesheetEle data doesn't contain the m strand
          //I will use the graph search in the same sheet id
          //formsheet function I build the graph information in onesheetEle data structure

          formsheet(onesheetEle);

          for( w=0;w<onesheetEle.size();w++)
          {
             onesheetEle[w].setchecked(false);
          }
          sheetid++;
          sheetset.clear();
          for(int y=0;y<onesheetEle.size();y++)
          {   
    //       cout<<"eleInstru[y].getcheckResult() "<<eleInstru[y].getcheckResult()<<" ";
    //       cout<<endl;
            // if(onesheetEle[y].getNeibor().size() != 0 && onesheetEle[y].getcheckResult()==false)
             if(onesheetEle[y].getcheckResult()==false)
             {
                search(onesheetEle, y,sheetset,sheetid);
                for(idindex=0;idindex<onesheetEle.size();idindex++)
                {
                   if(onesheetEle[idindex].getcheckResult() == false)
                   {
                      //onesheetEle[idindex].print();
                      sheetid++;
                      //cout<<"the sheetid is "<<sheetid<<endl;
                      break;
                   }
                }
             }
          }
           /*
           cout<<"the sheetset size is "<<sheetset.size()<<endl;
           for(w=0;w<sheetset.size();w++)
               sheetset[w].print();
           cout<<"the sheet id content is "<<endl;
           for(w=0;w<onesheetEle.size();w++)
              onesheetEle[w].print();
           */
 
           resetSheetid(eleInstru,onesheetEle,sheetset,sheetid);

           //this step I will assign the sheet id information to eleInstru data structure
           /*
           for( w=0;w<onesheetEle.size();w++)
              cout<<onesheetEle[w].getsheetId()<<"   "<<endl;
            exit(0);
           */
        }
        //above code I reset the eleInstru sheetid  information
        /*
        for(w=0;w<eleInstru.size();w++)
        {
           if(eleInstru[w].getElmentType() == "E")
           cout<<eleInstru[w].getsheetId()<<" ";
        }
        exit(0);
        */
     }
     else
     {
          fillsheetdefault(eleInstru);
          sheetset = sheetsetdefault;
     }
    fillsheetdefault(eleInstru);
    sheetset = sheetsetdefault;
    /*
    cout<<"the option is "<<option<<endl;
    cout<<"sheetset size is "<<sheetset.size()<<endl;
    for(i=0;i<sheetset.size();i++)
    {
         sheetset[i].print();
    }
    exit(0);
    */
    for(i=0 ;i<size ;i++)
    {
       intMatrix[i][i] = '*';
       v1 = eleInstru[i].getVector();
       size1 = eleInstru[i].getresidueGroup().size();
       elecoord = new point [size1]; 
       begin1 = eleInstru[i].getElmentType();
       chainB1 = eleInstru[i].getChainId();
       for(k=0 ; k<size1;k++)
       {
           elecoord[k].xcoord = eleInstru[i].getresidueGroup()[k].get_X();
           elecoord[k].ycoord = eleInstru[i].getresidueGroup()[k].get_Y();
           elecoord[k].zcoord = eleInstru[i].getresidueGroup()[k].get_Z();
       }
       for(j=i+1 ;j<size ;j++)
       {
           end1 = eleInstru[j].getElmentType();
           chainE1 = eleInstru[j].getChainId();
           p1 = eleInstru[i].getVector().initPoint;
           p2 = eleInstru[i].getVector().endPoint;
           p3 = eleInstru[j].getVector().initPoint;
           p4 = eleInstru[j].getVector().endPoint;
           v2 = eleInstru[j].getVector();
           overlap = overlapBetweenVector(p1 , p2 ,p3 ,p4);
           //cout<<"overlap is "<<overlap<<endl;
           size2 = eleInstru[j].getresidueGroup().size();
           elecoord2 = new point [size2];
           for(k=0 ; k<size2;k++)
           {
               elecoord2[k].xcoord = eleInstru[j].getresidueGroup()[k].get_X();
               elecoord2[k].ycoord = eleInstru[j].getresidueGroup()[k].get_Y();
               elecoord2[k].zcoord = eleInstru[j].getresidueGroup()[k].get_Z();
           }
           distance = distanceBetweenTwoElm(elecoord,size1, elecoord2, size2);
           //cout<<"the distacen is "<<distance<<endl;        
           //DISTANCE_DEFAULT    11.0 // 8 A this parament need to be adjusted
           //OVERLAP_DEFAULT     2.5  // 2 A, adjusted to real world
           //ANGLE_DEFAULT       85.0 
    //     cout<<"the distance is "<<distance<<" overlap  is "<<overlap;
           anglevalue = angle(v1.initPoint.xcoord,v1.initPoint.ycoord,v1.initPoint.zcoord,
                              v1.endPoint.xcoord,v1.endPoint.ycoord,v1.endPoint.zcoord,
                              v2.initPoint.xcoord,v2.initPoint.ycoord,v2.initPoint.zcoord,
                              v2.endPoint.xcoord,v2.endPoint.ycoord,v2.endPoint.zcoord
                             );
           int testing = 0;
/*
           fprintf(bug,"segment %4d%c -- %4d%c distance %8.3f overlap %8.3f angle %8.3f\n",i+1,begin1[0],j+1,end1[0],distance,overlap, anglevalue);
*/
           /*
            if(i==0 && j==2)
            {
                 cout<<"the distance is "<<distance<<" overlap is "<<overlap<<endl;
                 cout<<"(eleInstru[i].getElmentType() "<<eleInstru[0].getElmentType()<<endl;
                 cout<<"eleInstru[i].getElmentType() "<<eleInstru[2].getElmentType()<<endl;
                 exit(0);
            }
           */
            if(overlap > OVERLAP_DEFAULT && distance <= DISTANCE_DEFAULT )
            {
              if(eleInstru[i].getElmentType()=="E" && eleInstru[j].getElmentType()=="E")
              {
                 paralorNot = false;
                 /*
                 if(i==0 && j==2)
                 {
                     cout<<"this both E section "<<endl;
                     eleInstru[0].print();
                     exit(0);
                 }
                 */
                 if(h_bond_E(eleInstru[i] ,eleInstru[j]  , paralorNot)==true)
                 {
                      if(paralorNot == false)
                      {
                          intMatrix[i][j] =  't';
                      }
                      else
                          intMatrix[i][j] = 'c';
                 }
                 else
                 {
                    if( checkInsameSheet(i , j,sheetset) == true)
                    {
                         intMatrix[i][j] = '-';
                    }
                    else
                    { 
                      if(anglevalue > 180 - ANGLE_DEFAULT)
                         intMatrix[i][j] =  'v';
                      else if(anglevalue < ANGLE_DEFAULT)
                          intMatrix[i][j] = 'u';
                      else
                         intMatrix[i][j] = 'N';
                    }
                 
                 }
               }
              if(eleInstru[i].getElmentType()!="E" || eleInstru[j].getElmentType()!="E") 
               {
                       if(anglevalue > 180 - ANGLE_DEFAULT)
                       {
                          intMatrix[i][j] = 'v';
                       }
                       else if(anglevalue < ANGLE_DEFAULT)
                       {
                          intMatrix[i][j] = 'u';
                       }
                       else
                          intMatrix[i][j] = 'N';
                }
              }
           else
               { 
                     paralorNot = false;
                     if(h_bond_E(eleInstru[i] ,eleInstru[j]  , paralorNot)==true)
                     {
                      //cout<<"paralorNot is "<<paralorNot<<endl;
                         if(paralorNot == false)
                         {
                            intMatrix[i][j] =  't';
                         }
                         else
                            intMatrix[i][j] = 'c';
                      }
                      else
                         intMatrix[i][j] = '-';
                 }
               
           delete [] elecoord2;
       } 
       delete [] elecoord;
     } 
   //cout<<"the interaction matrix in this positin   "<<intMatrix[20][22]<<endl;
   fprintf(matrixptr,"%-32s",pdbid);
   fprintf(matrixptr,"%4d",eleInstru.size()); 
   for(i = 0; i<eleInstru.size(); i++)
   {
   	  pseudoBeginpos=eleInstru[i].getbeginId();
   	  pseudoEndpos=eleInstru[i].getendId();
   	  thiseleLength=eleInstru[i].getlength();
   	  trueBeginpos=rpseudoAtomindex[pseudoBeginpos];
   	  trueEndpos=rpseudoAtomindex[pseudoEndpos];
   	  cout<<pseudoBeginpos<<"--"<<pseudoEndpos<<":"<<trueBeginpos<<"--"<<trueEndpos<<endl;
   	  //matrixptr<<pseudoBeginpos<<"--"<<pseudoEndpos<<":"<<trueBeginpos<<"--"<<trueEndpos<<endl;
      fprintf(matrixptr,"%c%c%5s--%5s %4d %8.3f%8.3f%8.3f%8.3f%8.3f%8.3f",eleInstru[i].getElmentType()[0],
              eleInstru[i].getChainId()[0],trueBeginpos.c_str(),trueEndpos.c_str(),thiseleLength,
              eleInstru[i].getFirstP().xcoord,eleInstru[i].getFirstP().ycoord,
              eleInstru[i].getFirstP().zcoord, eleInstru[i].getLastP().xcoord,
              eleInstru[i].getLastP().ycoord,eleInstru[i].getLastP().zcoord);
   }
   fprintf(matrixptr,"\n");
   if(sheetset.size() == 0)
   {
      fprintf(matrixptr,"sheet %4d\n",sheetset.size());
   }
   else
   {
       for(i=0;i<sheetset.size();i++)
       {
          fprintf(matrixptr,"sheet %6d",i+1);
          fprintf(matrixptr,"%6d",sheetset[i].getSheet().size());
          for(j=0;j<sheetset[i].getSheet().size();j++)
          {
             fprintf(matrixptr,"%6d",sheetset[i].getSheet()[j]+1);
          }
          fprintf(matrixptr,"\n");
       }
   }
   for(i= 0; i<eleInstru.size(); i++)
   {
      for(j= i;j<eleInstru.size();j++)
      {
           fprintf(matrixptr,"%c",intMatrix[i][j]);
      }
   }
   fprintf(matrixptr,"\n");
   for(i=0;i<size;i++)
       delete [] intMatrix[i];
   delete [] intMatrix;
   cout<<"after delete "<<endl;
}
//this function is used to get the element a and b if they are H_bond
//and the parameter paralorNot is true they are parallel otherwose is
//false.
void control::printInterMatrix(vector<element> &eleInstru , char **a, FILE *bug,char *pdbid)
{   
    int i = 0;
    int j;
    char chainId;
    for(int y=0;y<eleInstru.size();y++)
    {
       chainId = eleInstru[y].getChainId()[0];
       fprintf(bug,"segment %3d range %4d%c--%4d%c\n",y+1,eleInstru[y].getbeginId(),chainId,eleInstru[y].getendId(),chainId);
    }
    int size = eleInstru.size();
    for(i=0;i<size;i++)
       fprintf(bug,"%3d",i+1);
    fprintf(bug,"\n");
    for(i=0;i<size;i++)
    {
       for(j=0;j<i;j++)
          fprintf(bug,"%3c",' ');
       for(int z=i;z<size;z++)
       {
          fprintf(bug,"%3c",a[i][z]);
       }
       fprintf(bug,"   %2d%c",i+1,eleInstru[i].getElmentType()[0]);
       fprintf(bug,"\n");
    }
    fclose(bug);
}
           
    
       
bool control::h_bond_E(element a , element b , bool& paralorNot)
{
    int i , j ,k;
    vector<Fresidue> a1 , b1;
    vector<h_bond> H_bond1;
    vector<h_bond> residueSet;
    vector<int> indexSet;
    h_bond temp;
    int number = 0;
    if(a.getbool() == true && b.getbool() == true)
    { 
       a1 = a.getresidueGroup();
       b1 = b.getresidueGroup();
       
       for(i=0 ;i<a1.size();i++)
       {
          H_bond1  = a1[i].getH_bond_residue();
          /*
          if(a1[i].getRsidNumber() == 30)
          {
             cout<<"a1[i].getH_bond_residue() size is "<<a1[i].getH_bond_residue().size()<<endl;
             a1[i].print();
          } 
          */
          for(j=0;j<H_bond1.size();j++)
          {
              residueSet.push_back(H_bond1[j]);
          }
       }
       for(j=0;j<residueSet.size();j++)
          for(k=0;k<b1.size();k++)
          {
             if(residueSet[j].getRnum() == b1[k].getRsidNumber()&&
                residueSet[j].getChaiName()== b.getChainId())
             {
                 indexSet.push_back(k);
                 number++;
             }
          }
 
       if(number>=2)
       {
           if(indexSet[0] > indexSet[indexSet.size()-1])
               paralorNot = false;
           else
               paralorNot = true;
           return true;
       }
       else
           return false;
    }
    else 
           return false;
}
       
void control::prinLogfile(vector<element> &eleInstru , char * filename)    
{
   int i,j;
   char logfileName[50];
   FILE *fvet;
   char nextname[20];
   char preName[20] = "P_PDB";
   strcpy(logfileName,filename);
   strcat(logfileName,".log");
   fvet = fopen(logfileName,"w");
   cout<<"the size is "<<eleInstru.size()<<endl;
   if(fvet == NULL)
   {
       cout<<"the "<<logfileName<<" can't open"<<endl;
       exit(0);
   }
   char cpfile[20];
   i = 3;
   cout<<"filename is "<<filename<<endl;
   while(filename[i] != '\0')
   {
     cpfile[i-3] = filename[i];
     i++;
   }
   cout<<"the i is "<<i<<endl;
   cpfile[i-3] = '\0';
   cout<<"cpfile is "<<cpfile<<endl;
   //exit(0);
   strcpy(nextname ,cpfile);
   strcat(cpfile,".pdb");
   i=0;
   while(nextname[i] != '\0')
   {
     if(nextname[i]>=97 && nextname[i]<=122)
      nextname[i] = nextname[i] - 32;
      i++;
   }
   strcat(preName , nextname);
   strcat(preName , "TXT");
   fprintf(fvet,"Get Molecule PDB User %s %s Heteroatom -Reference_Object\n",cpfile,preName);
   fprintf(fvet,"Display Molecule Only Atoms Trace %s*\n",preName);
   fprintf(fvet,"Center -World %s Center_of_mass\n",preName);
   fprintf(fvet,"Color Molecule Atoms %s Specified Specification 247,247,247\n",preName);
   for(j=0;j<eleInstru.size();j++)
   {
      if(eleInstru[j].getElmentType() == "E")
      {
           fprintf(fvet,"Color Molecule Atoms %s:%c%d-%c%d:CA Specified Specification yellow\n",
                    preName,eleInstru[j].getChainId()[0],eleInstru[j].getbeginId(),eleInstru[j].getChainId()[0],eleInstru[j].getendId());
           fprintf(fvet,"Color Molecule Atoms %s:%c%d Specified Specification orange\n",preName,eleInstru[j].getChainId()[0],eleInstru[j].getbeginId());
      }
      if(eleInstru[j].getElmentType() == "H")
      {
           fprintf(fvet,"Color Molecule Atoms %s:%c%d-%c%d:CA Specified Specification cyan\n",
                    preName, eleInstru[j].getChainId()[0],eleInstru[j].getbeginId(),eleInstru[j].getChainId()[0],eleInstru[j].getendId());
           fprintf(fvet,"Color Molecule Atoms %s:%c%d Specified Specification orange\n",preName,eleInstru[j].getChainId()[0],eleInstru[j].getbeginId());
      }     
   } 
  fprintf(fvet,"Color Molecule Atoms %s:%c%d Specified Specification green\n",preName,eleInstru
[0].getChainId()[0],eleInstru[0].getbeginId());
   fclose(fvet);
}
    
void control::cpmyvector(myvector &a , myvector &b)
{
    b.initPoint.xcoord = a.initPoint.xcoord;
    b.initPoint.ycoord = a.initPoint.ycoord;
    b.initPoint.zcoord = a.initPoint.zcoord;
    b.endPoint.xcoord = a.endPoint.xcoord;
    b.endPoint.ycoord = a.endPoint.ycoord;
    b.endPoint.zcoord = a.endPoint.zcoord;
}

void control::printmyvector(myvector &a)
{
   cout<<"the first point "<<a.initPoint.xcoord<<" "<<a.initPoint.ycoord<<" "<<a.initPoint.zcoord;
   cout<<endl;
   cout<<"the last  point "<<a.endPoint.xcoord <<" "<<a.endPoint.ycoord <<" "<<a.endPoint.zcoord;
   cout<<endl;
}

void control::printPoint(point &a)
{
   cout<<"point is :"<<endl;
   cout<<a.xcoord<<" "<<a.ycoord<<" "<<a.zcoord<<endl;
}
//this function I get the a graph data structure that form a E sheet.
void control::formsheet(vector<element> &eleInstru)
{
   int i ,j;
   vector<int> neior;
   int Neisize;
   bool paralorNot = false;
   int mcount = 0;

   for(i=0;i<eleInstru.size();i++)
   {
      if(eleInstru[i].getdeleteStr() == 'm')
        mcount++;
   }
   /*
   for(i=0;i<eleInstru.size();i++)
   {
      eleInstru[i].print();
   }
   for(i=0;i<eleInstru.size();i++)
   {
      if(eleInstru[i].getElmentType() == "E")
      cout<<eleInstru[i].getsheetId()<<endl;
   }
   */ 
   for(i=0 ;i<eleInstru.size();i++)
   {
        for(j=0;j<eleInstru.size();j++)
        {
          if(i !=j)
          {
              //if function h_bond_E return value is true, it means that there are H-bond between
              //these two elements
              if(h_bond_E(eleInstru[i] ,eleInstru[j]  , paralorNot)==true)
              {
                  //cout<<" "<<i<<" neibor "<<j;
                  eleInstru[i].addNeibor(j);
              }
           }
        }
    }
   /*
   for(i=0 ;i<eleInstru.size();i++)
   {
        neior = eleInstru[i].getNeibor();
        Neisize = neior.size();
        if(Neisize == 0)
          cout<<i<<" don't have link"<<endl;
        else
        {

           for(j=0 ; j<Neisize ;j++)
           {
              cout<<i<<" have neibor :"<<neior[j]<<" ";
           }
        }
        cout<<endl;
   }

   exit(0); 
   */
   
}
//this function is used to check how many E element in the same sheet.
void control::search(vector<element> &eleInstru , int i,vector<sheet> &sheetset,int sheetid)
{
    //0 means color white, 1 means color gray,2 means color black
    int k,j;
    int *d , *parent,*color;
    int size;
    int u;
    sheet shetemp;
    queue<int> que;
    vector<int> neighbor;
    size = eleInstru.size();
    d = new int [size];
    parent = new int[size];
    color = new int [size];
    for(k=0;k<eleInstru.size();k++) 
    {
     if(k != i)
     {
       eleInstru[k].setColor(0);
       color[k] = 0;
       d[k] = 10000;
       parent[k] = -1;
     }
    }
    eleInstru[i].setColor(1);
    color[i] = 1;
    eleInstru[i].setsheetId(sheetid);
    eleInstru[i].setchecked(true);
    d[i] = 0;
    parent[i] = -1;
    que.push(i);
    while(que.size() != 0)
    {
      u = que.front();
//      cout<<"u value is "<<u<<endl;
      neighbor.clear();
      neighbor = eleInstru[u].getNeibor();
      for(j=0 ; j<neighbor.size();j++)
      {
//         cout<<"neighbor[j] "<<neighbor[j]<<endl;
//         cout<<"color neighbor[j] "<<color[neighbor[j]]<<endl;
         if(color[neighbor[j]] == 0)
         {
             eleInstru[neighbor[j]].setColor(1);
             color[neighbor[j]] = 1;
             eleInstru[neighbor[j]].setchecked(true);
             eleInstru[neighbor[j]].setsheetId(sheetid);
             d[neighbor[j]] = d[i] + 1;
             que.push(neighbor[j]);
         }
       }
         shetemp.addelement(que.front());
         que.pop();
         color[u] = 2;
         eleInstru[u].setColor(2);
         eleInstru[u].setchecked(true); 
    }
   sheetset.push_back(shetemp);
   delete [] d;
   delete [] parent;
}  

      
       
#endif


