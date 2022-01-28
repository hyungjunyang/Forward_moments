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

// Block.cpp: implementation of the Block class.

#include "Block.h"

Block::Block(double v1_1, double v1_2, double v2_1, double v2_2, 
double dx, double dy, double dz, int i, int j, int k, Point pntLeft, Point pntCent)
: vx1(v1_1), vx2(v1_2), vy1(v2_1), vy2(v2_2),
  dx_(dx), dy_(dy), dz_(dz), i_(i), j_(j), k_(k), pnt_lb(pntLeft), pnt_ct(pntCent) {
  debug = true;

  dvx = (vx2 - vx1 ) / dx_;
  dvy = (vy2 - vy1 ) / dy_;
  wellType_     = false;
  wellType_Prod = false;
  wellType_Injc = false;
  boundType_Prod = false;
  boundType_Injc = false;
  xzero = 1e-10;
}

Block::Block(double dx, double dy, double dz, int i, int j, int k, 
             Point pntLeft, Point pntCent)
: dx_(dx), dy_(dy), dz_(dz), i_(i), j_(j), k_(k), 
  pnt_lb(pntLeft), pnt_ct(pntCent) {
  
  dvx = (vx2 - vx1 ) / dx_;
  dvy = (vy2 - vy1 ) / dy_;
  wellType_     = false;
  wellType_Prod = false;
  wellType_Injc = false;
  boundType_Prod = false;
  boundType_Injc = false; 
  xzero = 1e-10;
}

Block::~Block(){
}

void Block::setVelocity(double v1_1, double v1_2, double v2_1, double v2_2){
  vx1 = v1_1, vx2 = v1_2, vy1 = v2_1, vy2 = v2_2;
  dvx = (vx2 - vx1 ) / dx_;
  dvy = (vy2 - vy1 ) / dy_;
}

void Block::getV(Point& p, double& vx, double& vy){
   double x_len = p.getX() - pnt_lb.getX();
   double y_len = p.getY() - pnt_lb.getY();
   
   if(x_len < 0. || y_len < 0. ) {
      cout << "wrong block!! " <<endl;
           p.PrintXY(cout);
      pnt_lb.PrintXY(cout);
      exit(0);
   }
   vx = vx1 + dvx * x_len;
   vy = vy1 + dvy * y_len;
}

//pnt_lb is at the left-bottom corner
void Block::getMoveInfo(Point& p, double& travel_time, int & ib, int &jb){
   debug = false;
   debug = true;
   double x_len = p.getX() - pnt_lb.getX();
   double y_len = p.getY() - pnt_lb.getY();
   if(fabs (x_len) < 1.e-10) x_len = 0;
   if(fabs (y_len) < 1.e-10) y_len = 0;
   if(x_len < 0. || y_len < 0. ) {
       cout << "wrong block!! " <<endl;
       cout <<" x_len = " << x_len;
       cout <<" y_len = " << y_len << endl;
            p.PrintXY(cout);
       pnt_lb.PrintXY(cout);
       exit(0);
   }
// assert(x_len > 0 || y_len > 0);

   if( debug ) {
      cout<<"particle location("<<p.getX()<<','<<p.getY()<<')';
      cout<<"; block("<< i_<<','<<j_<<')';
      cout<<"; block corner("<<pnt_lb.getX()<<','<<pnt_lb.getY()<<')';
      cout<<"; x_len = "<<x_len<<", y_len = "<<y_len<<endl;
   }

   double vx = vx1 + dvx * x_len;
   double vy = vy1 + dvy * y_len;
   
   if(debug) {
      cout<<" vx1 = "<<vx1<< " vx = "<<vx<< " vx2 = " << vx2 <<endl;
      cout<<" vy1 = "<<vy1<< " vy = "<<vy<< " vy2 = " << vy2 <<endl;
      cout<<" dx = " << dx_ << " dy = " << dy_  << endl;
      cout<<" dvx= " << dvx << " dvy = " << dvy <<endl;
   }
   double dt[4] = {-1.e16,-1.e16,-1.e16,-1.e16};

   if(vx1 < 0 && vx < 0){ 
      if (fabs(dvx) < xzero) 
         dt[0] = (0. - x_len)/vx1;
      else 
         dt[0] = log(vx1/vx)/dvx;
   }

   if(vx2 > 0 && vx > 0){
      if (fabs(dvx) < xzero) 
          dt[1] = (dx_ - x_len)/vx2;
      else
          dt[1] = log(vx2/vx)/dvx;
   }

   if(vy1 < 0 && vy < 0){
      if (fabs(dvy) < xzero)
         dt[2] = (0. - y_len)/vy1;
      else
         dt[2] = log(vy1/vy)/dvy;
   }

   if(vy2 > 0 && vy > 0){
      if (fabs(dvy) < xzero) 
         dt[3] = (dy_ - y_len)/vy2;
      else
         dt[3] = log(vy2/vy)/dvy;
   }

   if(debug){
      for(int i=0; i< 4; ++i) cout<<dt[i]<<' '<<endl;
   }

   //   determine the exit side

   double delt = 1e20;
   int iexit = -1;
   if(        dt[0] > 0. && fabs((dt[0]-dt[2])/dt[0]) <= xzero ) {
       delt = dt[0];
       iexit = 4;
   } else if( dt[0] > 0. && fabs((dt[0]-dt[3])/dt[0]) <= xzero ) {
       delt = dt[0];
       iexit = 5;
   } else if( dt[1] > 0. && fabs((dt[1]-dt[2])/dt[1]) <= xzero ) {
       delt = dt[1];
       iexit = 6;
   } else if( dt[1] > 0. && fabs((dt[1]-dt[3])/dt[1]) <= xzero ) {
       delt = dt[1];
       iexit = 7;
   } else {
     for(int i = 0; i < 4; ++i){
	 //if(dt[i] > 0 && dt[i] < delt ){
         //if(dt[i] >= 0 && dt[i] < delt ){
	 if(dt[i] >= 0 && dt[i] < delt){
            delt = dt[i];
            iexit= i;
         }
	 if(fabs(dt[i]) < 0.0000001){
            delt = 0;
            iexit= i;
         }
     }
   }

   if(debug) cout<<"dt = "<<delt<<" iexit = "<<iexit<<endl;
   if(iexit<0){
     cout<<"ERROR"<<endl;
   }
   assert( iexit >= 0 );
    
// calculate the exit coordinate
   double  deltx, delty;

   if(fabs(dvx) > xzero ) 
      deltx = (vx * exp(dvx * delt)-vx1)/dvx;
   else 
      deltx = x_len + vx1 * delt;
   
   if(deltx > dx_) deltx = dx_ - 1e-10;

   if(fabs(dvy) > xzero )
      delty = (vy * exp(dvy*delt)-vy1)/dvy;
   else
      delty = y_len + vy1 * delt;

   if(delty > dy_) delty = dy_ - 1e-10;

// time of flying in this cell (cumulative time from the beginning)

   travel_time += delt;

// next block number
   switch( iexit ){
      case 0:
         deltx = 0.;
         ib = i_ - 1; 
         jb = j_;
         break;
      case 1:
         deltx = dx_;
         ib = i_ + 1;
         jb = j_;
         break;
      case 2:
         delty = 0.;
         ib = i_;
         jb = j_ - 1;
         break;
      case 3:
         delty = dy_;
         ib = i_;
         jb = j_ + 1;
         break;
      case 4:
	 deltx = 0;
	 delty = 0;
	 ib = i_ - 1;
         jb = j_ - 1;
	 break;
      case 5:
	 deltx = 0;
	 delty = dy_;
	 ib = i_ - 1;
         jb = j_ + 1;
	 break;
      case 6:
	 deltx = dx_;
	 delty = 0;
	 ib = i_ + 1;
         jb = j_ - 1;
	 break;
      case 7:
	 deltx = dx_;
	 delty = dy_;
	 ib = i_ + 1;
         jb = j_ + 1;
	 break;	 
      default:
         cout<<"Wrong at determining the exit face!"<<endl;
         exit(0); 
   }
//
   if(debug){
       cout<<"dx = " << deltx<<", dy = "<<delty<<endl;
       cout<<"Particle location for next block is: (";
          cout<<pnt_lb.getX()+deltx<<','<<pnt_lb.getY()+delty<<')'<<endl<<endl;
   }
     
   p = Point(pnt_lb.getX() + deltx, pnt_lb.getY() + delty,ib,jb);
}

void Block::changeVeloSign() { 
     vx1 = -vx1;
     vx2 = -vx2;
     vy1 = -vy1;
     vy2 = -vy2;
     dvx = -dvx;
     dvy = -dvy;
}

void Block::switchProdInjc(){
     bool tmp = wellType_Injc;
     wellType_Injc = wellType_Prod;
     wellType_Prod = tmp;
}

void Block::PrintBlock(ostream & os){
  os <<" i  = "<<i_  <<" j = " << j_ ;
  os <<" dx = "<<dx_ <<" dy = "<< dy_ ;
  os <<" vx1 = "<<vx1<<" vx2 = " << vx2 <<" vy1 = "<< vy1 <<" vy2 = "<<vy2
     <<" welltype = "<< wellType_Injc <<' '<<wellType_Prod<<endl;
  // pnt_lb.PrintXY(os);
}

