#ifndef HANDNESS_H
#define HANDNESS_H
class handness
{
   private :
      vector<int> elementId;
      char    hand;
   public  :
      handness();
      void addElmentId(int id);
      vector<int> & getIdgroup();
      void  sethandness(char a);
      char  gethandness();
      void  print();
      void  clear();
};

handness::handness()
{
    ;
}

void handness::addElmentId(int id)
{
   elementId.push_back(id);
}

void handness::print()
{
   int i;
   for(i=0;i<elementId.size();i++)
      cout<<elementId[i]<<" ";
   cout<<hand<<endl;
}

vector<int> & handness::getIdgroup()
{
   return elementId;
}
void handness::clear()
{
   elementId.clear();
   hand = ' ';
}
void handness::sethandness(char a)
{
   hand = a;
}

char handness::gethandness()
{
   return hand;
}
#endif
