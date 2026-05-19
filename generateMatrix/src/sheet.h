#ifndef SHEET_H
#define SHEET_H
class sheet{
   private :
         vector<int> eleid;
         int Nosheet;
   public :
          sheet();
          void addelement(int eleid);
          vector<int> &getSheet();
          void print();
          void setsheetid(int id);
          ~sheet();
};

void sheet::setsheetid(int id)
{
  Nosheet = id;
}

sheet::sheet()
{
  ;
}
sheet::~sheet()
{
   ;
}

void sheet::print()
{
   cout<<"this sheet contain the element"<<endl;
   int i=0;
   for(i=0 ; i<eleid.size();i++)
     cout<<eleid[i]<<" ";
   cout<<endl;
   ;
}

void sheet::addelement(int eleNo)
{
   eleid.push_back(eleNo);
}

vector<int> & sheet::getSheet()
{
   return eleid;
}

#endif


        
