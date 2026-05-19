#ifndef ELEMENT_H
#define ELEMENT_H 
#include "external.h"
class element
{
 private:
   vector<Fresidue> residueGroup;
   string chainId;
   string elementType;
   myvector elementVector;
   vector<point> project;
   int whichelement;
   bool havebond;
   vector<int>  neiborElent;
   bool checked;
   int  color;
   int  length;
   int  beginNum;
   int  endNum;
   int  sheetid;
   char deleteStr ; // m is short one strand that should be deleted, space is the strand that should be kept
 public :
   element();
   element(string cName , string elmentType);
   ~element();
   void setsheetId(int id);
   int  getsheetId();
   bool operator<(const element &a) const;
   bool operator==(const element &a) const;
   vector<Fresidue> & getresidueGroup();
   string getChainId();
   string  getElmentType();
   myvector getElementVector();
   void addResidue(Fresidue &a);
   void setwhichelement(int number);
   void setVector(myvector &a);
   int getwhichElement();
   void print();
   int  getResidueSize();
   void setEltype(string a);
   myvector &getVector();
   bool  getbool();
   void addProPoint(point a);
   vector<point> & getProject(); 
   void addNeibor(int n);
   vector<int> & getNeibor();
   void setchecked(bool a);
   bool  getcheckResult();
   void  setColor(int k);
   int   getColor();
   int  getbeginId();
   int  getendId();
   point getFirstP();
   point getLastP();   
   void setStarEndnum(int start,int end);
   void setlength(int lengh);
   int  getlength();
   int  getbeginNum();
   int  getendNum();
   int  getopresidue();
   void setdeleteStr(char letter);
   char getdeleteStr();
   bool withelement(int residueNum, char chainId);
};

bool element::withelement(int renum, char chid)
{
   bool find = false;

   char insert = ' ';
   char insertion1 = ' ';
   char insertion2 = ' ';
   //cout<<"the renum is "<<renum<<endl;
   //cout<<"the chid  is "<<chid<<endl;
   //cout<<"the insert is "<<insert<<endl;
   if(renum > beginNum && chid == chainId[0])
   {
     if(renum < endNum)
     {
        find = true;
     }
     if(renum == endNum)
     {
        if(insert <= insertion2)
          find = true;
        else
          find = false;
     }
     if(renum > endNum)
       find = false;
   }
   //this case I conder the inserstion case,the first condition
   //means that this residue in behind the first residue
   if(renum == beginNum && chid ==chainId[0] && insert >= insertion1)
   {
      if(renum < endNum && chid == chainId[0])
         find =  true;
      if(renum == endNum && chid == chainId[0])
      {
         if(insert <= insertion2)
            find = true;
         else
            find = false;
      }
      if(renum > endNum && chid == chainId[0])
         find = false;
   }
   if(renum == endNum && chid == chainId[0] && insert <= insertion2)
     find = true;
   //cout<<"the find is "<<find<<endl;
   return find;
}
void element::setdeleteStr(char letter)
{
   deleteStr = letter;
}

char element::getdeleteStr()
{
   return deleteStr;
}

void element::setsheetId(int id)
{
    sheetid = id;
}

int  element::getsheetId()
{
    return sheetid;
}

int element::getbeginNum()
{
   return beginNum;
}

int element::getopresidue()
{
   if(residueGroup.size()>0)
   {
      return residueGroup[residueGroup.size()-1].getRsidNumber();
   }
}
int element::getendNum()
{
   return endNum;
}
point element::getFirstP()
{
   point a;
   a = project[0];
   return a;
}

bool element::operator==(const element &a) const
{
  if(chainId == a.chainId && beginNum == a.beginNum && endNum == a.endNum)
     return true;
  else
     return false;
}

bool element::operator<(const element &a) const
{
   if(chainId > a.chainId)
      return false;
   else if(chainId == a.chainId)
   {
      if(beginNum > a.beginNum)
        return false;
      else
        return true;
   }
   else
      return true;
}

void element::setStarEndnum(int start,int end)
{
   beginNum = start;
   endNum = end;
}

void element::setlength(int lengh)
{
   length = lengh;
}

int element::getlength()
{
   return length;
}
point element::getLastP()
{
   point a;
   int size;
   size = project.size();
   a = project[size-1];
   return a;
}
int element::getbeginId()
{
   int id;
   id = residueGroup[0].getRsidNumber();
   return id;
}

int element::getendId()
{
   int id;
   int size;
   size = residueGroup.size();
   id = residueGroup[size-1].getRsidNumber();
}
void element::setchecked(bool a)
{
    checked = a;
}
void element::setColor(int k)
{
   color = k;
}

int element::getColor()
{
    return color;
}
bool element::getcheckResult()
{
    return checked;
}
void element::addNeibor(int n)
{
   neiborElent.push_back(n);
}

vector<int> & element::getNeibor()
{
     return neiborElent;
}
void element::setEltype(string a)
{
   elementType = a;
}

void element::addProPoint(point a)
{
   project.push_back(a);
}

vector<point> & element::getProject()
{
   return project;
}
bool element::getbool()
{
    return havebond;
}
void element::setVector(myvector &a)
{
   
    elementVector.initPoint.xcoord = a.initPoint.xcoord;
    elementVector.initPoint.ycoord = a.initPoint.ycoord;
    elementVector.initPoint.zcoord = a.initPoint.zcoord;
    elementVector.endPoint.xcoord = a.endPoint.xcoord;
    elementVector.endPoint.ycoord = a.endPoint.ycoord;
    elementVector.endPoint.zcoord = a.endPoint.zcoord;
}
element::~element()
{
  ;
}

element::element()
{
  havebond = false ;
  checked = false;
  sheetid = 0;
  deleteStr = ' ';
}

myvector& element::getVector()
{
/*
   elementVector.initPoint = project[0];
   elementVector.endPoint = project[project.size()-1];
*/
   return elementVector;
}
int element::getResidueSize()
{
    return residueGroup.size();
}
void element::print()
{
   int j=0;
   /*
   for(int i=0;i<residueGroup.size();i++)
   {
     cout<<"the residue "<<residueGroup[i].getRsidNumber()<<" have ";
     for(j=0;j<residueGroup[i].getH_bond_residue().size();j++)
     {
         cout<<"the residue "<<residueGroup[i].getRsidNumber()<<" have ";
         cout<<"the chain is "<<residueGroup[i].getH_bond_residue()[j].getChaiName()<<" ";
         cout<<"the number is "<<residueGroup[i].getH_bond_residue()[j].getRnum()<<" ";
         cout<<endl;
     } 
   }
   cout<<"the project size is "<<project.size()<<endl;
   cout<<endl;
   */
   cout<<"the size is "<<residueGroup.size()<<" "; 
   cout<<"sheetid is "<<sheetid<<" ";
   cout<<"the beginNum is "<<beginNum<<"  "<<"end number is "<<endNum<<" chain name is "<<chainId<<endl;
}

element::element(string cName , string eleType)
{
    chainId = cName;
    elementType = eleType;
    deleteStr = ' ';
}

vector<Fresidue> & element::getresidueGroup()
{
    return residueGroup;
}

string element::getChainId()
{
   return chainId;
}

string element::getElmentType()
{
   return elementType;
}

void element::addResidue(Fresidue &a)
{
    if(a.getbool() == true)
      havebond  = true;
    residueGroup.push_back(a);
}

void element::setwhichelement(int number)
{
    whichelement = number;
}

int element::getwhichElement()
{
    return whichelement;
}



#endif
   
