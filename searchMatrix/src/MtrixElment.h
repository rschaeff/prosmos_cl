#ifndef MTRIXELEMENT_H
#define MTRIXELEMENT_H
#include "external.h"
#include <string.h>
class matrixElment
{
    private :
     char elementType;
     char chainName;
     point initpoint;
     point endpoint;
     string  beginNum;
     int  lenght;
     string  endNum;
    public :
        matrixElment();
        void setElmentType(char tpe);
        void setChainName(char name);
        void setFpoint(point b);
        void setLpoint(point a);
        char getEltype();
        char getChainName();
        point getFpoint();
        point getLpoint();
        void print();
        void setlength(int len);
        void setBeginNum(string begin);
        void setEndNum(string endNum);
        string  getBeginNum();
        string  getEndNum();
        int  getlenght();
};
void matrixElment::setlength(int len)
{
   lenght = len;
}
matrixElment::matrixElment()
{
     ;
}

char matrixElment::getEltype()
{
   return elementType;
}
char matrixElment::getChainName()
{
   return chainName;
}
int matrixElment::getlenght()
{
   //return endNum - beginNum+1;//need my change
   return lenght;
}

void matrixElment::setBeginNum(string begin)
{
    beginNum = begin;
}

void matrixElment::setEndNum(string endNum1)
{
    endNum = endNum1;
}

string matrixElment::getBeginNum()
{
    return beginNum;
}

string matrixElment::getEndNum()
{
   return endNum;
} 
void matrixElment::print()
{
    cout<<"the element type is "<<elementType<<" the chain name is "<<chainName
        <<"the begin number is "<<beginNum<<" the end number    is "<<endNum<<endl;
    cout<<"the initial point of vector "<<endl<<" xcoord is "<<initpoint.xcoord<<" ycoord "
        <<initpoint.ycoord<<" zcoord is "<<initpoint.zcoord<<endl;
    cout<<"the last point of vector "<<endl<<" xcoord is "<<endpoint.xcoord<<" ycoord "
        <<endpoint.ycoord<<" zcoord is "<<endpoint.zcoord<<endl<<endl<<endl;
}
void matrixElment::setElmentType(char tpe)
{
    elementType = tpe;
}

point matrixElment::getFpoint()
{
    return initpoint;
}


point matrixElment::getLpoint()
{
    return endpoint;
}

void matrixElment::setChainName(char name)
{
    chainName = name;
}

void matrixElment::setFpoint(point  b)
{
    initpoint.xcoord = b.xcoord;
    initpoint.ycoord = b.ycoord;
    initpoint.zcoord = b.zcoord;
}
void matrixElment::setLpoint(point  b)
{
    endpoint.xcoord = b.xcoord;
    endpoint.ycoord = b.ycoord;
    endpoint.zcoord = b.zcoord;
}
#endif

