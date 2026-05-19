#ifndef FPASS_H
#define FPASS_H
class fpass{
   private :
       vector<elecol> parent;
       int lenth;
       elecol  currentNode;
   public :
       fpass();
       fpass(elecol a);
       vector<elecol> &getparent();
       int    getLenth();
       elecol  getcurrentNode();
       void  addNodeToparent(elecol a);
       void  inheritPass(fpass a);
       void  print();
};
void fpass::inheritPass(fpass a)
{
    int i;
    for(i=0;i<a.parent.size();i++)
    {
       parent.push_back(a.parent[i]);
    }
}

fpass::fpass()
{
   lenth = 0;
}

void fpass::print()
{
   int i;
   cout<<"the current node is "<<endl;
   currentNode.print();
   for(i=0;i<parent.size();i++)
   {
      cout<<"this is parent node "<<i<<endl;
      parent[i].print();
   }
} 
void fpass::addNodeToparent(elecol a)
{
    parent.push_back(a);
    lenth++;
}
fpass::fpass(elecol a)
{
   currentNode = a;
}
vector<elecol>& fpass::getparent()
{
   return parent;
}

int fpass::getLenth()
{
   return lenth;
}

elecol fpass::getcurrentNode()
{
  return currentNode;
}
#endif       
