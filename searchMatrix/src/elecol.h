#ifndef ELECOL_H
#define ELECOL_H
#include <string>
using namespace std;
class elecol{
  private:
     char eleTyp;
     int col;
     char chainame;
     //this parameter is used to store neighbor information for the E in the same sheet
     vector<int> neighbor;
     char paraOranti;
     string color;
     int parent;
  public :
     elecol();
     elecol(char a , int col1,char chainN);
     char geteleTyp();
     char getchain();
     void setparalleloranti(char a);
     char getparalleloranti();
     int    getcol();
     void   seteleTyp(char a);
     void   setAllvalue(char a,int col1);
     void   setcol(int col);
     void   print();
     vector<int>  getneighbor();
     void addneighbor(int a);
     void setcolorparent(string color , int par);
     string getcolor();
     void setcolor(string a);     
     int getparent();
     char getparal();
};
char elecol::getparal()
{
   return paraOranti;
}
int elecol::getparent()
{
   return parent;
}

void elecol::setcolor(string a)
{
   color = a;
}

string elecol::getcolor()
{
    return color;
}

void elecol::setcolorparent(string colr , int par)
{
    color = colr;
    parent = par;
}
vector<int>  elecol::getneighbor()
{
   return neighbor;
}

void elecol::setparalleloranti(char a)
{
   paraOranti = a;
}

char elecol::getparalleloranti()
{
   return paraOranti;
}    
void elecol::addneighbor(int a)
{
   neighbor.push_back(a);
}
char elecol::getchain()
{
    return chainame;
}
void elecol::setAllvalue(char a, int col1)
{
    eleTyp = a;
    col = col1;
} 
elecol::elecol()
{
    ;
}

void elecol::print()
{
  // if(neighbor.size()>0)
  {
   int i;
   cout<<"the element type is "<<eleTyp
   <<"the col  number  is "<<col<<" the chain is "<<chainame<<endl;
   cout<<"its neighbor is "<<endl;
   
   for(i=0 ; i<neighbor.size();i++)
   {
      cout<<neighbor[i]<<" ";
   }
   cout<<endl;
  
   cout<<"the color is "<<color<<" parent is "<<parent<<endl;
  }
}

elecol::elecol(char a , int col1,char chainN)
{
   eleTyp = a;
   col = col1;
   chainame = chainN;
}

char elecol::geteleTyp()
{
   return eleTyp;
}

int  elecol::getcol()
{
   return col;
}
void elecol::seteleTyp(char a)
{
   eleTyp = a;
}

void  elecol::setcol(int col1)
{
   col = col1;
} 
#endif
