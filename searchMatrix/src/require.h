//this file is used to category the requirement type when the user 
//submit their different requirements
//there are several different requirements
//1 sheetD 1 4 that means 1 and 4 element is in the different sheet
//2 sheetS 1 2 3 4 that means that 1 2 3 4 elements are in the same sheet
//3 chainD 1 4 that means 1 and 4 element is in the different chains
//4 chainS 1 ,2 3,4 means that 1,2,3,4 elements are in the same chain
//5 length 1 E 4 - 8 means that element 1 B-srand its length is between 4 and 8 residues.
#ifndef REQUIRE_H
#define REQUIRE_H
class require
{
    private :
          string requireType;
          vector<int> elementSet;
          //these four attributes is used for the length require
          int  eleId;
          int  startRange;
          int  endRange;
          string eleType;
    public :
          require();
          void setRequireType(string type);
          void setEleType(string type);
          void setStartEndRange(int start, int end);
          void setEleId(int id);
          void addElement(int id);
          vector<int> &getElementSet();
          string getrequireType();
          //the following foure function is used for the length requirement
          int geteleId();
          int getstartRange();
          int getendRange();
          string geteleType();
          void  print();
};
require::require()
{
   ;
}
void require::print()
{
    int i;
    if(requireType == "length")
    {
       cout<<"the element id is "<<eleId<<"the element type is"<<eleType<<" the start range is "<<startRange<<" the end range is "<<endRange<<endl;
    }
    else
    {
       cout<<"the require type is "<<requireType<<endl;
       cout<<"the element set is ";
       for(i=0;i<elementSet.size();i++)
        cout<<elementSet[i]<<" ";
    }
    cout<<endl; 
}
string require::getrequireType()
{
   return requireType;
}

void require::setRequireType(string type)
{
    requireType = type;
}

void require::setEleType(string type)
{
    eleType = type;
}

void require::setStartEndRange(int start, int end)
{
    startRange = start;
    endRange = end;
}

void require::setEleId(int id)
{
   eleId = id;
}

void require:: addElement(int id)
{
   elementSet.push_back(id);
}

vector<int> & require::getElementSet()
{
   return elementSet;
}   
int require::geteleId()
{
   return eleId;
}

int require::getstartRange()
{
    return startRange;
}

int require::getendRange()
{
    return endRange;
}
string require::geteleType()
{
     return eleType;
}
#endif

 
