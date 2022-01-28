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
/*           --       wellArc.h       --
             --    class definition    --

Note: here assume that wellare within the block

blockPtr : a block pointer, to which wellarc points

len_     : the number of points used to discretize the wellarc;
ptSet_   : the set of points used to discretize the wellarc;
theta    : the angles at each point
weight   : intgration weights
vx, vy   : vx and vy at each point
v_total  : total velocity at each point
radius_  : the radius of the well arc;

*/

#ifndef WELLARCH
#define WELLARCH

#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include "Point.h"
#include "Block.h"

const double Pi = 3.1415926535897932;
void gauleg(double x1, double x2, double x[], double w[], int n);

class WellArc{
public:
   WellArc(Block**, Block*, bool, bool, int, int, int, int, int, int,
           double, double, double, double, double);
   WellArc(const WellArc & a) { copy(a); }       // copy constructor
   
  ~WellArc(){free();}
   void allocMem ();

   void free();
   WellArc & operator=(const WellArc &);   // assignment
   Point& operator[](int);
   const Point& operator[](int) const;

   int length() const { return len_; }
   bool isProdWell()   { return wellType_Prod;}
   bool isInjcWell()   { return wellType_Injc;}

   double  getRadius() {return radius_ ;}
   double* getTheta()  {return theta ;  }
   double* getWeight() {return weight ; }
   double* getTotalV() {return v_total ;}
   int getI() {return i_;}
   int getJ() {return j_;}
   int getK() {return k_;}
   int getId() {return id;}

   // set part
   void setLength(int len) {len_ = len;}
   void setPointOrder(int nx, int ny);
   
   void switchProdInjc();

   //Launching Type
   void CalcCirclePoint(int nx, int ny, double x0, double y0);
   void FaceType(int nx, int ny, int i0, int j0);
   void setPoint(int, int, int, double xtemp, double ytemp);

   void PrintArc(ostream & os);
   void DisplayArc();
   Block* getBlock() { return wellBlockPtr;}
private:
   int i_, j_, k_, id;

   Block** blockPtr;
   Block* wellBlockPtr;
   
   double smallValue;
   int  len_;                         // length of WellArc
   Point* ptSet_;
   double *theta, *weight, *vx, *vy, *v_total;
   double radius_;
   bool wellType_Injc;
   bool wellType_Prod;
   void copy(const WellArc &);
};

#endif
