//////////////////////////////////////////////////
//                                              //
//                 Liyong Li                    //
//      Reservoir Simulation Research Team      //
//        Chevron Petroleum Technology Co.      //
//          1300 Beach Blvd., Rm. 3166          //
//            La Habra, CA 90631-6374           //
//             Phone:(562) 694-7366             //
//              Fax:(562) 694-7565              //
//            Email liyl@chevron.com            //
//                                              //
//                 Version 1.                   //
//                Apr. 07, 1999                 //
//////////////////////////////////////////////////

//             --       Block.h          --
//             --    class definition    --

#ifndef BLOCKH
#define BLOCKH

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include "Common.h"
#include "Point.h"

class Block{
public:
   Block(double, double, double, double, 
         double dx, double dy, double dz, int i, int j, int k, Point, Point);
   Block(double, double, double, int, int, int, Point, Point);
   virtual ~Block();
   int    getI()   {return pnt_lb.getI();}
   int    getJ()   {return pnt_lb.getJ();}
   int    getK()   {return pnt_lb.getK();}
   int    getNode(){return node_;}
   double getDx()  {return dx_;}
   double getDy()  {return dy_;}
   double getDz()  {return dz_;}
   double getVx1() {return vx1;}
   double getVx2() {return vx2;}
   double getVy1() {return vy1;}
   double getVy2() {return vy2;}
   Point getBlockCtPnt() {return pnt_ct;} 
   Point getBlockLbPnt() {return pnt_lb;}
   void setProdWell()  { wellType_Prod = true;}
   void setInjcWell()  { wellType_Injc = true;}
   void setWellBlock() { wellType_ = true;}

   void setInjcBound()  { boundType_Injc = true;} 
   void setProdBound()  { boundType_Prod = true;} 

   bool isProdWell()   { return wellType_Prod;}
   bool isInjcWell()   { return wellType_Injc;}
 
   bool isProdBound()   { return boundType_Prod;}
   bool isInjcBound()   { return boundType_Injc;}

   bool wellBlock()    { return wellType_;}
   void switchProdInjc();
   void changeVeloSign();
   void getMoveInfo(Point &, double &, int &, int &);

   void setVelocity(double, double, double, double);
   void getV(Point& p, double&, double&);
   
   void PrintBlock(ostream & os);

private:           
   bool debug;

   Point pnt_lb;     //point of left bottom in the cell
   Point pnt_ct;     //point of center in the cell
   double vx1;
   double vx2;
   double vy1;
   double vy2;
   double dvx;
   double dvy;
   double dx_;
   double dy_;
   double dz_;
   double xzero;
   int i_;
   int j_;
   int k_;
   int node_;
   bool wellType_;
   bool wellType_Injc;
   bool wellType_Prod;

   bool boundType_Prod; 
   bool boundType_Injc; 

};

#endif
