#ifndef H_BOND_H 
#define H_BOND_H
#include <iostream>
using namespace std;
class h_bond
{
   private :
             int masterNum;
             string masterChain;
             int residueNum;
             string chainName;
   public  :
             h_bond();
             h_bond(string Cname , int rnum);
             void setRnum(int Rnum);
             void setchainName(char cName);
             int getRnum();
             string getChaiName();
             void print();
             void setMaster(string cname,int rnum);
             void setSuborder(string cname,int rum);
             bool operator==(const h_bond &) const;
             int  getMasterNum();
             string getMasterChain();
             bool operator<=(const h_bond &a) const;
             bool operator>(const h_bond &a) const;
             const h_bond &operator=(const h_bond &a) ;
};
const h_bond& h_bond::operator=(const h_bond &a)
{
  masterNum = a.masterNum;
  masterChain = a.masterChain;
  residueNum = a.residueNum;
  chainName = a.chainName;
}
bool h_bond::operator==(const h_bond &a) const
{
    if(residueNum == a.residueNum && chainName == a.chainName&&masterNum == a.masterNum && masterChain==a.masterChain)
       return true;
    else
       return false;
}
               
bool h_bond::operator>(const h_bond &a) const
{
   if(masterChain == a.masterChain)
   {
      if(masterNum>a.masterNum)
      {
          return true;
      }
      else
          return false;
   }
   if(masterChain<a.masterChain)
      return false;
   if(masterChain > a.masterChain)
      return true;
} 
bool h_bond::operator<=(const h_bond &a) const
{
   if(masterChain == a.masterChain)
   {
      if(masterNum<=a.masterNum)
      {
          return true;
      }
      else
          return false;
   }
   if(masterChain<a.masterChain)
      return true;
   if(masterChain > a.masterChain)
      return false;
}
  




void h_bond::setSuborder(string cname,int rnum)
{
   chainName = cname;
   residueNum = rnum;
} 
int h_bond::getMasterNum()
{
    return masterNum;
}

string h_bond::getMasterChain()
{
    return masterChain;
}
void h_bond::setMaster(string cname,int rnum)
{
    masterNum = rnum;
    masterChain = cname;
}

h_bond::h_bond()
{
   ;
}
h_bond::h_bond(string Cname , int rnum)
{
       chainName = Cname;
       residueNum = rnum;
}
      
void h_bond::print()
{
   cout<<"the master num is "<<masterNum;
   cout<<"    the master chain  "<<masterChain;
   cout<<"    the resid num  is   "<<residueNum;
   cout<<"    the chain name is "<<chainName<<endl;
}

void h_bond::setRnum(int Rnum)
{
    residueNum = Rnum;
}

void h_bond::setchainName(char cName)
{
   chainName = cName;
}

int h_bond::getRnum()
{
    return residueNum;
}

string h_bond::getChaiName()
{
    return chainName;
}

#endif

     
