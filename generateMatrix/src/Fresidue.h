#ifndef FRESIDUE_H
#define FRESIDUE_H
#include "h_bond.h"
#include <string>
using namespace std;
class Fresidue
{  
 public :
	   // default constuctor
	       Fresidue();
	       Fresidue(const Fresidue &r);
	       ~Fresidue();
	       string &getRsidName();
               Fresidue(string cname ,int Rnum,string Rname ,double x ,double y , double z, string eleNam);
	       int getRsidNumber();
	       void printResidue();
	       string getchaiName();
	       void setchaiName(string chaiName);
	       void setResidName(string name);
               void setredueNum(int redNum);
	       const Fresidue &operator=(const Fresidue &);
               void set_X(double X);
               void set_Y(double Y);
               void set_Z(double Z);
               double get_X();
               double get_Y();
               double get_Z();
	       void setbool(bool a);
               bool getbool();	   
               void print();
               void addHbondResidue(h_bond &a);
               void setEletyp(string eletype);
               vector<h_bond>& getH_bond_residue();
  protected :
	        string residueName;
	     	int residueNumber;
		string chainName;
                double X_coordinate;
                double Y_coordinate;
                double Z_coordinate; 
               // double distanItoI3;
               // double TorSionAngle;
                string elementType;
                vector<h_bond> h_bond_residue;
                bool havebond;
};


/*******************************************************************************/
double Fresidue::get_X()
{
   return X_coordinate;
}

double Fresidue::get_Y()
{
   return Y_coordinate;
}

double Fresidue::get_Z()
{
   return Z_coordinate;
}

vector<h_bond>&  Fresidue::getH_bond_residue()
{
    return h_bond_residue;
}
Fresidue::Fresidue()
{
    havebond = false  ;
}

void Fresidue::setEletyp(string eletype)
{
    elementType = eletype;
}

Fresidue::Fresidue(string cname ,int Rnum,string Rname ,double x ,double y , double z, string eleNam)
{
    
   residueName = Rname;
   residueNumber = Rnum;
   chainName = cname;
   X_coordinate = x;
   Y_coordinate = y;
   Z_coordinate = z; 
   //distanItoI3  = disiToi3;
   //TorSionAngle = Torsion;
   elementType  = eleNam;
}
Fresidue::Fresidue(const Fresidue & r)
{
    
   residueName = r.residueName;
   residueNumber = r.residueNumber;
   chainName = r.chainName;
   X_coordinate = r.X_coordinate;
   Y_coordinate = r.Y_coordinate;
   Z_coordinate = r.Z_coordinate; 
   //distanItoI3  = r.distanItoI3;
   //TorSionAngle = r.TorSionAngle;
   elementType  = r.elementType;
   h_bond_residue = r.h_bond_residue; 
}

const Fresidue & Fresidue::operator=(const Fresidue &r)
{   
   residueName = r.residueName;
   residueNumber = r.residueNumber;
   chainName = r.chainName;
   X_coordinate = r.X_coordinate;
   Y_coordinate = r.Y_coordinate;
   Z_coordinate = r.Z_coordinate;
   //distanItoI3  = r.distanItoI3;
   //TorSionAngle = r.TorSionAngle;
   elementType  = r.elementType;
   h_bond_residue = r.h_bond_residue;
}
    
void Fresidue::print()
{
  /*
   cout<<"the chain name     is "<<chainName<<endl;
   cout<<"the residue number is "<<residueNumber<<endl;
   cout<<"the residue name   is "<<residueName<<endl;
   cout<<"the coordinate     is "<<X_coordinate<<"  "<<Y_coordinate<<"  "<<Z_coordinate<<endl;
   cout<<"the element type   is "<<elementType<<endl;
  */
   int i;
   for(i=0;i<h_bond_residue.size();i++)
      h_bond_residue[i].print();
   cout<<endl;
   
        
}


void Fresidue::addHbondResidue(h_bond &a)
{
    int i=0;
    bool judge = false;
    for(i=0;i<h_bond_residue.size();i++)
    {
       if(h_bond_residue[i].getRnum()==a.getRnum()&&h_bond_residue[i].getChaiName()==a.getChaiName())
       {
          judge = true; 
       }
    }
    if(judge == false)
    {
          h_bond_residue.push_back(a);
          havebond = true;
    }
}
Fresidue::~Fresidue()
{
	;
}

int Fresidue::getRsidNumber()
{
	return residueNumber;
}


	


string Fresidue::getchaiName()
{
	return chainName;
}

string & Fresidue::getRsidName()
{
	return residueName;
}

void Fresidue::setchaiName(string chaiName)
{
	chainName = chaiName;
}

void Fresidue::setResidName(string name)
{
	residueName = name;
}
void Fresidue::setredueNum(int redNum)
{
	residueNumber = redNum;
}
void Fresidue::setbool(bool a)
{
     havebond = a;
}
bool Fresidue::getbool()
{
	return havebond;
}

#endif


	

