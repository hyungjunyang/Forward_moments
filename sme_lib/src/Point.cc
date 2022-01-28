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

#include "Point.h"

Point& Point::operator = (const Point & rhs){
  if(this == &rhs) return *this;
   x_ = rhs.x_;
   y_ = rhs.y_;
   z_ = rhs.z_;
   i_ = rhs.i_;
   j_ = rhs.j_;
   k_ = rhs.k_;

  return *this;
}

void Point::PrintXY(ostream & os){
  os.precision(9);
  os.setf(ios::fixed, ios::floatfield);
  /* 
  os << setw(12) << x_ << " " << setw(12) << y_ <<endl;
  */
  os << " i = " << i_ << " j= " << j_ ;
  os << " x = " << x_ << " y= " << y_ << endl; 
  /*
  os << " i = "<<i_<<" j= " << j_ <<" k= "<<k_;
  os << " x = "<<x_<<" y= " << y_ <<" z= "<<z_ <<endl; 
  */
}

