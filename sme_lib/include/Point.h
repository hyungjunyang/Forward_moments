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

//             --       Point.h          --
//             --    class definition    --

#ifndef POINTH
#define POINTH

#include <iostream>
#include <math.h>
#include <assert.h>
#include <fstream>
#include <iomanip>

using namespace std;

class Point{


public:
   Point():x_(0.), y_(0.), z_(0),i_(0), j_(0),k_(0) {}        
   Point(double x, double y, int i, int j) 
   : x_(x), y_(y), i_(i), j_(j){z_=0.;k_=0;}
   Point(double x, double y, double z) : x_(x), y_(y), z_(z){i_=j_=k_=0;}
   Point(double x, double y, double z, int i, int j, int k) 
	 : x_(x), y_(y), z_(z), i_(i), j_(j), k_(k){}
   Point(const Point &rhs)
	 : x_(rhs.x_),y_(rhs.y_), z_(rhs.z_),i_(rhs.i_),j_(rhs.j_), k_(rhs.k_){}
   virtual ~Point(){ ;}

//OVERLOADED OPERATORS:
   Point& operator= (const Point &rhs);  
   double getX() {return x_;}
   double getY() {return y_;}
   double getZ() {return z_;}
   int getI() {return i_;}
   int getJ() {return j_;}
   int getK() {return k_;}
   
   void PrintXY(ostream & os);	

   private:    
   double x_, y_, z_;
   int    i_, j_, k_;
};

#endif
