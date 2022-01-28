/////////////////////////////////////////////////
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

#include "Particle.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Particle::Particle(Point p, double time_, Domain& dmn ):
 particle_point(p),
 particle_time_(0),
 domain_(dmn)
{
   try_Cartesian = true;
   //try_Cartesian = false;

   particle_block = domain_.point2block(p);
   
   // p.PrintXY(cout);
   // particle_block->PrintBlock(cout);
   
   trajectory.push_back(particle_point);
   time_flight.push_back(particle_time_);
}

Particle::~Particle(){
   delete [] node;
   delete [] node_i;
   delete [] node_j;
   delete [] costh;
   delete [] sinth;
   delete [] ds;
   delete [] x_x1;
   delete [] y_y1;
   delete [] x2_x;
   delete [] y2_y;

   delete [] dv_deta;
   delete [] v_avg;
   delete [] v_var;

   delete [] tau_Avg0;
   delete [] tau_Avg2;
   delete [] tau_Vari;
   
   delete [] tau_Avg0Trans;
   delete [] tau_Avg2Trans;
   delete [] tau_VariTrans;

   delete [] tau_AvgPdf;

   //cartesian
   if(try_Cartesian) {
      delete [] tau_Avg0Carte;
      delete [] tau_Avg2Carte;
      delete [] tau_VariCarte;
   }
}

void Particle::initialization() {
   node   = new int[length];
   node_i = new int[length];
   node_j = new int[length];
   costh  = new double[length];
   sinth  = new double[length];
   ds     = new double[length]; 

   x_x1   = new double[length];
   y_y1   = new double[length];
   x2_x   = new double[length];
   y2_y   = new double[length];
   dv_deta  = new double[length];
   v_avg    = new double[length];
   v_var    = new double[length];

   tau_Avg0     = new double[length];
   tau_Avg2     = new double[length];
   tau_Vari     = new double[length];
   
   // Transformed system
   tau_Avg0Trans = new double[length];
   tau_Avg2Trans = new double[length];
   tau_VariTrans = new double[length];
   
   // pdf  
   tau_AvgPdf   = new double[length];
   
   // Cartesian
   if(try_Cartesian) {
      tau_Avg0Carte = new double[length];
      tau_Avg2Carte = new double[length];
      tau_VariCarte = new double[length];
   }
}

Block* Particle::getBlock(){
   return particle_block;
}

void Particle::setBlock(Block* b){
   particle_block = b;
}

void Particle::move(){
   int i = 0, j = 0;
   particle_block->getMoveInfo(particle_point, particle_time_,i,j);

 //  particle_point.PrintXY(cout);

   particle_block = domain_.getBlock(i, j);
   trajectory.push_back( particle_point );
   time_flight.push_back(particle_time_);
}

/* Function: moveBackward
 * Author: Pipat Likanapaisal
 * --------------------------
 * The function track a particle along negative velocity field.
 * Be careful that the particle_point is at upstream, also
 * trajectory and time_flight are in backtracking fashion.
 */
void Particle::moveBackward()
{
	int i = 0, j = 0;
	particle_block->changeVeloSign();
	particle_block->getMoveInfo(particle_point, particle_time_,i,j);
	particle_block->changeVeloSign();
	particle_block = domain_.getBlock(i, j);
	trajectory.push_back(particle_point);
	time_flight.push_back(particle_time_);
}

/* Function: reorderData
 * Author: Pipat Likanapaisal
 * --------------------------
 * The function reorder trajectory and time_flight to forward
 * tracking fashion. This function is needed for backward
 * tracking method so that it can use other functions.
 * However, the current position is not changed.
 */
void Particle::reorderData()
{
	int vsize = trajectory.size();
	vector<Point> temptraj(trajectory);
	vector<double> temptime(time_flight);
	trajectory.clear();
	time_flight.clear();
	double sumtime = 0;
	double ptime = temptime.back();
	for(int ii=0;ii<vsize;ii++)
	{
	  //if(ii!=0){cout<<temptime.back()<<" "<<ptime<<" "<<(temptime.back()==ptime)<<endl;}
		if((ii!=0)&&(temptime.back()==ptime))
		{
			temptraj.pop_back();
			temptime.pop_back();
		}
		else
		{
			trajectory.push_back(temptraj.back());
			temptraj.pop_back();
			time_flight.push_back(sumtime+(ptime-temptime.back()));
			sumtime = time_flight.back();
			ptime = temptime.back();
			temptime.pop_back();
		}
	}
}

void Particle::CoordTrans() {
   length = trajectory.size();
   initialization();
   Block *block_i;
   double dx, dy, x1, y1, x2, y2;
   double x_p, y_p, x_po, y_po;
   double x_len, y_len, x_len2, y_len2;
   double vx, vx2, vx1, dvx; 
   double vy, vy2, vy1, dvy;

   int node_i1, node_j1, node_ij; 
   double var11_x , var11_x1,   var11_x2,   var11_x12;
   double var22_y , var22_y1,   var22_y2,   var22_y12;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double *Cv1v1   = domain_.getCv1v1();
   double *Cv2v2   = domain_.getCv2v2();
   double *Cv1v2   = domain_.getCv1v2();
   int num_blocks  = domain_.getNumBlocks();

   for(int i = 0; i < length; ++i) {
       block_i  = domain_.point2block( trajectory[i] );
       dx       = block_i->getDx();
       dy       = block_i->getDy();
       x1       = block_i->getBlockLbPnt().getX();
       y1       = block_i->getBlockLbPnt().getY();
       x2       = x1 + dx;
       y2       = y1 + dy;
       x_p      = trajectory[i].getX();
       y_p      = trajectory[i].getY();
       x_len    = x_p - x1;
       y_len    = y_p - y1;
       x_len2   = x2 - x_p;
       y_len2   = y2 - y_p;
       x_x1[i] = x_len  / dx;
       y_y1[i] = y_len  / dy;
       x2_x[i] = x_len2 / dx;
       y2_y[i] = y_len2 / dy;
       
       if(i == 0) {
          ds[0] = 0.;
       } else {  
          x_po = trajectory[i - 1].getX();
          y_po = trajectory[i - 1].getY();
          ds[i] = sqrt(  (x_p - x_po)*(x_p - x_po) 
                       + (y_p - y_po)*(y_p - y_po) 
                      );
       }

       vx1      = block_i->getVx1();
       vx2      = block_i->getVx2();
       vy1      = block_i->getVy1();
       vy2      = block_i->getVy2();
       dvx      = (vx2 - vx1 ) / dx;
       dvy      = (vy2 - vy1 ) / dy;
       vx       = vx1  + dvx * x_len;
       vy       = vy1  + dvy * y_len;
       v_avg[i] = sqrt(vx * vx + vy * vy); 
       costh[i]  = vx / v_avg[i];
       sinth[i]  = vy / v_avg[i];
       
       node_i1  = domain_.getIndex( block_i->getBlockLbPnt().getI() - 1,
                                    block_i->getBlockLbPnt().getJ()
                                  );
       node_j1  = domain_.getIndex( block_i->getBlockLbPnt().getI(),
                                    block_i->getBlockLbPnt().getJ() - 1
                                  );
       node_ij  = domain_.getIndex( block_i->getBlockLbPnt().getI(),
                                    block_i->getBlockLbPnt().getJ()
                                  );
       node[i]   = node_ij;
       node_i[i] = node_i1;
       node_j[i] = node_j1;
   
       //Pipat
	   int ii = block_i->getI();
	   int jj = block_i->getJ();
	   if((ii==0)||(jj==0)||(ii==domain_.getNx()-1)||(jj==domain_.getNy()-1))
	   {
		   var11_x1   = 0;
		   var11_x2   = 0;
           var11_x12  = 0;
           var22_y1   = 0;
           var22_y2   = 0;
           var22_y12  = 0;
           var12_x1y1 = 0;
           var12_x2y1 = 0;
           var12_x1y2 = 0;
           var12_x2y2 = 0;
	   }
	   else
	   {
		   var11_x1   = Cv1v1[ node_i1 + node_i1 * num_blocks ];
		   var11_x2   = Cv1v1[ node_ij + node_ij * num_blocks ];
           var11_x12  = Cv1v1[ node_i1 + node_ij * num_blocks ];
           var22_y1   = Cv2v2[ node_j1 + node_j1 * num_blocks ];
           var22_y2   = Cv2v2[ node_ij + node_ij * num_blocks ];
           var22_y12  = Cv2v2[ node_j1 + node_ij * num_blocks ];
           var12_x1y1 = Cv1v2[ node_i1 + node_j1 * num_blocks ];
           var12_x2y1 = Cv1v2[ node_ij + node_j1 * num_blocks ];
           var12_x1y2 = Cv1v2[ node_i1 + node_ij * num_blocks ];
           var12_x2y2 = Cv1v2[ node_ij + node_ij * num_blocks ];
	   }
	   //Pipat
	   /*Pipat
	   var11_x1   = Cv1v1[ node_i1 + node_i1 * num_blocks ];
       var11_x2   = Cv1v1[ node_ij + node_ij * num_blocks ];
       var11_x12  = Cv1v1[ node_i1 + node_ij * num_blocks ];
       var22_y1   = Cv2v2[ node_j1 + node_j1 * num_blocks ];
       var22_y2   = Cv2v2[ node_ij + node_ij * num_blocks ];
       var22_y12  = Cv2v2[ node_j1 + node_ij * num_blocks ];
       var12_x1y1 = Cv1v2[ node_i1 + node_j1 * num_blocks ];
       var12_x2y1 = Cv1v2[ node_ij + node_j1 * num_blocks ];
       var12_x1y2 = Cv1v2[ node_i1 + node_ij * num_blocks ];
       var12_x2y2 = Cv1v2[ node_ij + node_ij * num_blocks ];
	   Pipat*/
       
       if(block_i->getBlockLbPnt().getI() == 0 )
          var11_x1 = var11_x12 = var12_x1y1 = var12_x1y2 = 0.;
       if(block_i->getBlockLbPnt().getJ() == 0 ) 
          var22_y1 = var22_y12 = var12_x2y1 = var12_x1y1 = 0.;

       var11_x = var11_x2   * x_x1[i] * x_x1[i]
               + var11_x1   * x2_x[i] * x2_x[i]
               + var11_x12  * x_x1[i] * x2_x[i] * 2.;
       var22_y = var22_y2   * y_y1[i] * y_y1[i] 
               + var22_y1   * y2_y[i] * y2_y[i] 
               + var22_y12  * y_y1[i] * y2_y[i] * 2.;
       var12_xy= var12_x2y2 * x_x1[i] * y_y1[i]
               + var12_x1y1 * x2_x[i] * y2_y[i]
               + var12_x1y2 * x2_x[i] * y_y1[i]
               + var12_x2y1 * x_x1[i] * y2_y[i];
       
       v_var[i] = var11_x  * costh[i] * costh[i] 
             +    var22_y  * sinth[i] * sinth[i] 
             +2.* var12_xy * sinth[i] * costh[i];
       dv_deta[i] = (-dvx + dvy) * sinth[i] * costh[i];
   }
}

void Particle::CoordTrans2(double *Vxi_Vxi, double *Vxi_Eta, 
		           double *Eta_Eta){
   double *Veta_Veta= new double[length * length];
   double *Vxi_Veta = new double[length * length];

   double *Cv1v1   = domain_.getCv1v1();
   double *Cv2v2   = domain_.getCv2v2();
   double *Cv1v2   = domain_.getCv1v2();
   
   double var11_x , var11_x1,   var11_x2,   var11_x12,  var11_x21;
   double var22_y , var22_y1,   var22_y2,   var22_y12,  var22_y21;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double var21_xy, var21_x1y1, var21_x1y2, var21_x2y1, var21_x2y2;
   double C_vxvx, C_vyvy, C_vxvy, C_vyvx;

   int num_blocks  = domain_.getNumBlocks();
   Block *block_l1, *block_l2;

   int l1l2;   
   for(int l2 = 0; l2 < trajectory.size(); ++l2) {
       block_l2  = domain_.point2block( trajectory[l2] );
       //trajectory[l2].PrintXY(cout);

       for(int l1 = 0; l1 < trajectory.size(); ++l1) {
           block_l1  = domain_.point2block( trajectory[l1] );

           var11_x2   = Cv1v1[ node[l1]   + node[l2]   * num_blocks ];
           var22_y2   = Cv2v2[ node[l1]   + node[l2]   * num_blocks ];
           var12_x2y2 = Cv1v2[ node[l1]   + node[l2]   * num_blocks ];
           var21_x2y2 = Cv1v2[ node[l2]   + node[l1]   * num_blocks ];

           if(block_l1->getBlockLbPnt().getI() == 0 ) {
              var11_x12  = var12_x1y2 = 0;
           } else {
              var11_x12  = Cv1v1[ node_i[l1] + node[l2]   * num_blocks ];
              var12_x1y2 = Cv1v2[ node_i[l1] + node[l2]   * num_blocks ];
           }
             
           if(block_l2->getBlockLbPnt().getI() == 0 ) {
              var11_x21  = var21_x2y1 = 0;
           } else {
              var11_x21  = Cv1v1[ node[l1]   + node_i[l2] * num_blocks ];
              var21_x2y1 = Cv1v2[ node_i[l2] + node[l1]   * num_blocks ];
           }

           if(block_l1->getBlockLbPnt().getI() == 0 || 
              block_l2->getBlockLbPnt().getI() == 0 ) {
              var11_x1   = 0;
           } else {
              var11_x1   = Cv1v1[ node_i[l1] + node_i[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getJ() == 0 ) {
              var22_y12  = var21_x1y2 = 0;           
           } else {
              var22_y12  = Cv2v2[ node_j[l1] + node[l2]   * num_blocks ];
              var21_x1y2 = Cv1v2[ node[l2]   + node_j[l1] * num_blocks ];
           }
           
           if(block_l2->getBlockLbPnt().getJ() == 0 ) {
              var22_y21  = var12_x2y1 = 0;           
           } else {
              var22_y21  = Cv2v2[ node[l1]   + node_j[l2] * num_blocks ];
              var12_x2y1 = Cv1v2[ node[l1]   + node_j[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getJ() == 0 ||
              block_l2->getBlockLbPnt().getJ() == 0 ) {
              var22_y1   = 0;
           } else {
              var22_y1   = Cv2v2[ node_j[l1] + node_j[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getI() == 0 ||
              block_l2->getBlockLbPnt().getJ() == 0 ) {
              var12_x1y1 = 0.;  
           } else {
              var12_x1y1 = Cv1v2[ node_i[l1] + node_j[l2] * num_blocks ];
           }
           if(block_l1->getBlockLbPnt().getJ() == 0 ||
              block_l2->getBlockLbPnt().getI() == 0 ) {
              var21_x1y1 = 0.;
           }else {
              var21_x1y1 = Cv1v2[ node_i[l2] + node_j[l1] * num_blocks ];
           }
           
           C_vxvx = var11_x2   * x_x1[l1] * x_x1[l2]
                  + var11_x12  * x2_x[l1] * x_x1[l2]
                  + var11_x21  * x_x1[l1] * x2_x[l2]
                  + var11_x1   * x2_x[l1] * x2_x[l2];
           C_vyvy = var22_y2   * y_y1[l1] * y_y1[l2]
                  + var22_y12  * y2_y[l1] * y_y1[l2]
                  + var22_y21  * y_y1[l1] * y2_y[l2]
                  + var22_y1   * y2_y[l1] * y2_y[l2];
           C_vxvy = var12_x2y2 * x_x1[l1] * y_y1[l2]
                  + var12_x1y2 * x2_x[l1] * y_y1[l2] 
                  + var12_x2y1 * x_x1[l1] * y2_y[l2]
                  + var12_x1y1 * x2_x[l1] * y2_y[l2];
           C_vyvx = var21_x2y2 * y_y1[l1] * x_x1[l2]
                  + var21_x1y2 * y2_y[l1] * x_x1[l2]
                  + var21_x2y1 * y_y1[l1] * x2_x[l2]
                  + var21_x1y1 * y2_y[l1] * x2_x[l2];
           l1l2   = l1 + l2 * length;
           Vxi_Vxi[ l1l2 ]  = C_vxvx * costh[l1] * costh[l2]
                            + C_vyvx * sinth[l1] * costh[l2]
                            + C_vxvy * costh[l1] * sinth[l2]
                            + C_vyvy * sinth[l1] * sinth[l2];
           Veta_Veta[l1l2 ] = C_vxvx * sinth[l1] * sinth[l2]
                            - C_vyvx * costh[l1] * sinth[l2]
                            - C_vxvy * sinth[l1] * costh[l2]
                            + C_vyvy * costh[l1] * costh[l2];
           Vxi_Veta[ l1l2 ] =-C_vxvx * costh[l1] * sinth[l2]
                            - C_vyvx * sinth[l1] * sinth[l2]
                            + C_vxvy * costh[l1] * costh[l2]
                            + C_vyvy * sinth[l1] * costh[l2];
	   //cout<<l1<<" "<<l2<<" "<<C_vxvx<<" "<<C_vyvy<<" "<<C_vxvy<<" "<<C_vyvx<<endl;
       }
   }
   int l1l2m1;
   for(int l1 = 0; l1 < length; ++l1) {
       l1l2 = l1;
       Vxi_Eta[l1l2] = 0.;
   }

   for(int l2 = 1; l2 < length; ++l2) {
       for(int l1 = 0; l1 < length; ++l1) {
          l1l2   = l1 + l2       * length;
          l1l2m1 = l1 + (l2 - 1) * length;
          Vxi_Eta[l1l2] = Vxi_Eta[l1l2m1] 
                        +(   Vxi_Veta[l1l2]/v_avg[l2] 
                           + Vxi_Veta[l1l2m1]/v_avg[l2-1]
                         ) /2. * ds[l2];
       }
   }

   int l1m1l2, l1m1l2m1;
   for(int l2 = 0; l2 < length; ++l2) {
       l1l2 = l2 * length;
       Eta_Eta[l1l2] = 0.;
   }
   for(int l1 = 0; l1 < length; ++l1) {
       l1l2     = l1;
       Eta_Eta[l1l2] = 0.;
   }
   double tmp;
   for(int l2 = 1; l2 < length; ++l2) {
       for(int l1 = 1; l1 < length; ++l1) {
          l1l2     = l1     + l2 * length;   
          l1m1l2   = l1 - 1 + l2 * length;
          l1l2m1   = l1     + (l2 - 1) * length;   
          l1m1l2m1 = l1 - 1 + (l2 - 1) * length;
          tmp = ( Veta_Veta[l1l2]    /v_avg[l1  ]/v_avg[l2  ] 
                + Veta_Veta[l1m1l2]  /v_avg[l1-1]/v_avg[l2  ]
                + Veta_Veta[l1l2m1]  /v_avg[l1  ]/v_avg[l2-1]
                + Veta_Veta[l1m1l2m1]/v_avg[l1-1]/v_avg[l2-1]
                )/4.*ds[l1]*ds[l2];
          Eta_Eta[l1l2] = Eta_Eta[l1m1l2] + Eta_Eta[l1l2m1] 
                        + tmp -  Eta_Eta[l1m1l2m1];
       }
   }
   //==== checking ===
   /*
   int l1l1, l2l2;
   for(int l2 = 0; l2 < length; ++l2) {
       for(int l1 = 0; l1 < length; ++l1) {
           l1l2     = l1 + l2 * length;
	   l1l1     = l1 + l1 * length;
	   l2l2     = l2 + l2 * length;
	   double stdv = sqrt(Eta_Eta[l1l1] * Eta_Eta[l2l2] );
	   if(stdv > 0.00001) {
	      cout << l1 <<' ' << l2<<' ' << Eta_Eta[l1l2]/stdv << endl; 
	   }
       }
       cout<<endl;
   }
   exit(0);
   */
   //=== end of check ===
   delete [] Vxi_Veta;
   delete [] Veta_Veta;
}

void Particle::calcEtaDU(double *etadu) {

   int num_blocks   = domain_.getNumBlocks();
   double *dxCv1v1  = domain_.getDxCv1v1();  
   double *dyCv1v1  = domain_.getDyCv1v1();  
   double *dxCv2v2  = domain_.getDxCv2v2();  
   double *dyCv2v2  = domain_.getDyCv2v2(); 
   double *dx1Cv1v2 = domain_.getDx1Cv1v2();
   double *dx2Cv1v2 = domain_.getDx2Cv1v2();
   double *dy1Cv1v2 = domain_.getDy1Cv1v2(); 
   double *dy2Cv1v2 = domain_.getDy2Cv1v2();

   double tm1, tm2, sum;
   for(int i = 0; i < trajectory.size(); ++i) {
       sum = 0.;
       for(int j = 1; j <= i; ++j) {
           int k = j - 1;
	   /*
           tm1  = sinth[i] * costh[i] * sinth[j] * dxCv1v1 [ node[i] + node[k] * num_blocks ];
           tm1 += sinth[i] * sinth[i] * sinth[j] * 0.;
           tm1 -= costh[i] * costh[i] * sinth[j] * 0.;
           tm1 -= sinth[i] * costh[i] * sinth[j] * dy2Cv1v2[ node[k] + node[i] * num_blocks ];
           tm1 -= sinth[i] * costh[i] * costh[j] * dx1Cv1v2[ node[i] + node[k] * num_blocks ];
           tm1 -= sinth[i] * sinth[i] * costh[j] * 0.;
           tm1 += costh[i] * costh[i] * costh[j] * 0.;
           tm1 += sinth[i] * costh[i] * costh[j] * dyCv2v2 [ node[i] + node[k] * num_blocks ];
	   */
           tm1  = sinth[i] * costh[i] * sinth[k] * dxCv1v1 [ node[i] + node[k] * num_blocks ];
           tm1 += sinth[i] * sinth[i] * sinth[k] * 0.;
           tm1 -= costh[i] * costh[i] * sinth[k] * 0.;
           tm1 -= sinth[i] * costh[i] * sinth[k] * dy2Cv1v2[ node[k] + node[i] * num_blocks ];
           tm1 -= sinth[i] * costh[i] * costh[k] * dx1Cv1v2[ node[i] + node[k] * num_blocks ];
           tm1 -= sinth[i] * sinth[i] * costh[k] * 0.;
           tm1 += costh[i] * costh[i] * costh[k] * 0.;
           tm1 += sinth[i] * costh[i] * costh[k] * dyCv2v2 [ node[i] + node[k] * num_blocks ];

           tm1 /= v_avg[k];
           
           tm2  = sinth[i] * costh[i] * sinth[j] * dxCv1v1 [ node[i] + node[j] * num_blocks ];
           tm2 += sinth[i] * sinth[i] * sinth[j] * 0.;
           tm2 -= costh[i] * costh[i] * sinth[j] * 0.;
           tm2 -= sinth[i] * costh[i] * sinth[j] * dy2Cv1v2[ node[j] + node[i] * num_blocks ];
           tm2 -= sinth[i] * costh[i] * costh[j] * dx1Cv1v2[ node[i] + node[j] * num_blocks ];
           tm2 -= sinth[i] * sinth[i] * costh[j] * 0.;
           tm2 += costh[i] * costh[i] * costh[j] * 0.;
           tm2 += sinth[i] * costh[i] * costh[j] * dyCv2v2 [ node[i] + node[j] * num_blocks ];
	   
           tm2 /= v_avg[j];

           sum += (tm1 + tm2)/2.0 * ds[j];
       }
       etadu[i] = sum/ v_avg[i];
   }
}

// == TavAvg ==
void Particle::calcTauAvg(){

   calcTauAvg0();
   calcTauAvg2();
   //calcTauAvg2Trans();
   //calcTauAvgPdf();
   
   if(try_Cartesian) {
      try_Carte();
   }
}

// ===== Cartesian System =====
void Particle::try_Carte() {
   
   double * c_eta_eta = new double[length * length];
   double * c_eta_vx  = new double[length * length];
   double * kernel    = new double[length * length];
   double * c_eta_dvx = new double[length * length];
   double * c_vx_vx   = new double[length * length];
   double * c_vy_vx   = new double[length * length];
   double * dvx_deta  = new double[length ];
   double * delta_x   = new double[length ];
   double * delta_y   = new double[length ];
   double * vx_avg    = new double[length ];
   double * vy_avg    = new double[length ];

   calcEtaEtaCarte( c_eta_eta, kernel, c_eta_vx, c_eta_dvx, c_vx_vx,  c_vy_vx, dvx_deta,
		    delta_x, delta_y, vx_avg, vy_avg);
   calcTauAvg0Carte(vx_avg, delta_x);
   calcTauAvg2Carte(c_eta_eta, c_eta_vx, c_eta_dvx, dvx_deta, c_vx_vx, delta_x, delta_y, vx_avg);
   calcTauVariCarte(vx_avg, c_vx_vx, c_eta_eta, dvx_deta, c_eta_vx, delta_x);
                                
   delete [] c_eta_eta;
   delete [] c_eta_vx;
   delete [] c_eta_dvx;
   delete [] kernel;
   delete [] c_vx_vx;
   delete [] c_vy_vx;
   delete [] dvx_deta;
   delete [] delta_x;
   delete [] delta_y;
   delete [] vx_avg;
   delete [] vy_avg;
}

void Particle::calcTauAvg0Carte(double* vx, double* delta_x){
   
   double xzero = 1e-10;
   double vx_tmp, vx2, vx1, dvx, dx, x_len, dtx = 0.; 
   Block *block_i, *block_j;
   double t_sum = 0., dt = 0.;
   /*
   for(int i = 0; i < trajectory.size(); ++i) {
       block_i = domain_.point2block( trajectory[i] );
       x_len = trajectory[i].getX() 
             - block_i->getBlockLbPnt().getX();
       dx    = block_i->getDx();
       vx1   = block_i->getVx1();
       vx2   = block_i->getVx2();
       
       dvx   = (vx2 - vx1 ) / dx;
       vx_tmp   = vx1 + dvx * x_len;

       t_sum += dt;
       tau_Avg0Carte[i] = t_sum;
       
       if(vx1 < 0 && vx_tmp < 0){ 
          if (fabs(dvx) < xzero) 
             dtx = (0. - x_len)/vx1;
          else 
             dtx = log(vx1/vx_tmp)/dvx;
       }
       if(vx2 > 0 && vx_tmp > 0){
         if (fabs(dvx) < xzero) 
             dtx = (dx - x_len)/vx2;
         else
             dtx = log(vx2/vx_tmp)/dvx;
       }
       
       if(i + 1 < trajectory.size() ) {
          dt = dtx;
       }
   }
   */
   tau_Avg0Carte[0] = 0;
   for(int i = 1; i < trajectory.size(); ++i) {
       block_i = domain_.point2block( trajectory[i] );
       x_len = trajectory[i].getX() 
             - block_i->getBlockLbPnt().getX();
       dx    = block_i->getDx();
       vx1   = block_i->getVx1();
       vx2   = block_i->getVx2();

       dvx   = (vx2 - vx1 ) / dx;
       vx_tmp   = vx1 + dvx * x_len;

       if(vx1 < 0 && vx_tmp < 0){ 
          if (fabs(dvx) < xzero) 
             dtx = (0. - x_len)/vx1;
          else 
             dtx = log(vx1/vx_tmp)/dvx;
       }
       if(vx2 > 0 && vx_tmp > 0){
         if (fabs(dvx) < xzero) 
             dtx = (dx - x_len)/vx2;
         else
             dtx = log(vx2/vx_tmp)/dvx;
       }
       tau_Avg0Carte[i] = tau_Avg0Carte[i - 1] + dtx;
   }
   
   /*
   tau_Avg0Carte[0] = 0;
   for(int i = 1; i < trajectory.size(); ++i) {
       block_i = domain_.point2block( trajectory[i] );
       x_len = trajectory[i].getX() 
             - block_i->getBlockLbPnt().getX();
       dx    = block_i->getDx();
       vx1   = block_i->getVx1();
       vx2   = block_i->getVx2();
       dvx   = (vx2 - vx1 ) / dx;
       dt  = log(fabs(vx[i]/vx[i-1]))/fabs(dvx);
       tau_Avg0Carte[i] = tau_Avg0Carte[i - 1] + dt;
   }
   */
   /*
   tau_Avg0Carte[0] = 0;
   double tmp, tmp1, tmp2;
   for(int i = 1; i < trajectory.size(); ++i) {
       tmp1 = 1. / vx[i  ];
       tmp2 = 1. / vx[i-1];
       tmp  = (tmp1 + tmp2) * delta_x[i] /2.;	   
       tau_Avg0Carte[i] = tau_Avg0Carte[i - 1] + tmp;
   }*/
}

void Particle::calcEtaEtaCarte(double* c_eta_eta, double* kernel, double* c_eta_vx, double* c_eta_dvx,
                               double* c_vx_vx,  double* c_vy_vx, double* dvx_deta,
			       double* delta_x, double* delta_y, double* vx, double* vy) {

   Block *block_i, *block_l1, *block_l2;
	
   double dx, dy;
   double x1, y1, x2, y2;
   double x_p, x_po, y_p, y_po;
   double x_len, x_len2, y_len, y_len2;
   double vx1, vx2, vy1, vy2, dvx, dvy;
   
   // calculate   d <vx (x, eta )>    d <vx (x, eta )>     dx
   //           ---------------  = ----------------- -----------
   //               d <eta>              d x            d <eta>

   for(int i = 0; i < length; ++i) {
       block_i  = domain_.point2block( trajectory[i] );
       dx       = block_i->getDx();
       dy       = block_i->getDy();
       x1       = block_i->getBlockLbPnt().getX();
       y1       = block_i->getBlockLbPnt().getY();

       x_p      = trajectory[i].getX();
       y_p      = trajectory[i].getY();
       x_len    = x_p - x1;
       y_len    = y_p - y1;

       if(i == 0) {
          delta_x[0] = 0.;
	  delta_y[0] = 0.;
       } else {  
          x_po = trajectory[i - 1].getX();
	  y_po = trajectory[i - 1].getY();
          delta_x[i] = sqrt(  (x_p - x_po)*(x_p - x_po) );
	  delta_y[i] = sqrt(  (y_p - y_po)*(y_p - y_po) );
       }

       vx1      = block_i->getVx1();
       vx2      = block_i->getVx2();
       vy1      = block_i->getVy1();
       vy2      = block_i->getVy2();
       dvx      = (vx2 - vx1 ) / dx;
       dvy      = (vy2 - vy1 ) / dy;
       vx[i]    = vx1  + dvx * x_len;
       //cout << block_i->getBlockLbPnt().getI() << ' '  
       //     << block_i->getBlockLbPnt().getJ() << ' '
       //     << vx1 << ' ' << vx[i] << ' ' << vx2 << endl;
       vy[i]    = vy1  + dvy * y_len;
       dvx_deta[i] = dvx * vx[i]/vy[i];
       //cout << dvx_deta[i] << endl;
   }
   //cout << endl;
   //exit(0);
   //cout << endl;
   /*             
   for(int i = 0; i < length; ++i) {
       if(i == 0) {
          if( fabs(delta_y[i]) < 0.0000001) {
              dvx_deta[i] = 0.;		  
          } else {		  
	      dvx_deta[i] = (vx[i + 1] - vx[i])/delta_y[i];
	  }
       } else if( i == length - 1) {
          if( fabs(delta_y[i]) < 0.0000001) {
              dvx_deta[i] = 0.;		  
          } else {		  
	      dvx_deta[i] = (vx[i] - vx[i - 1])/delta_y[i];
	  }
       } else {
          if( fabs(delta_y[i]) < 0.0000001) {
              dvx_deta[i] = 0.;		  
          } else {	
	      dvx_deta[i] = (vx[i + 1] - vx[i - 1])/(delta_y[i]+delta_y[i-1]);	  
	  }
       }
       cout << dvx_deta[i] << endl;
   }
   */
   
   int iflag_algo = 1;
   double *Cv1v1   = domain_.getCv1v1();
   double *Cv2v2   = domain_.getCv2v2();
   double *Cv1v2   = domain_.getCv1v2();
   
   double var11_x , var11_x1,   var11_x2,   var11_x12,  var11_x21;
   double var22_y , var22_y1,   var22_y2,   var22_y12,  var22_y21;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double var21_xy, var21_x1y1, var21_x1y2, var21_x2y1, var21_x2y2;
   double C_vxvx, C_vyvy, C_vxvy, C_vyvx;

   int num_blocks  = domain_.getNumBlocks();

   int l1l2;   
   double f_old, f_new, f_old_dvx, f_new_dvx;
   for(int l2 = 0; l2 < trajectory.size(); ++l2) {
       block_l2  = domain_.point2block( trajectory[l2] );
       f_old     = 0.;
       for(int l1 = 0; l1 < trajectory.size(); ++l1) {
           block_l1  = domain_.point2block( trajectory[l1] );

           var11_x2   = Cv1v1[ node[l1]   + node[l2]   * num_blocks ];
           var22_y2   = Cv2v2[ node[l1]   + node[l2]   * num_blocks ];
           var12_x2y2 = Cv1v2[ node[l1]   + node[l2]   * num_blocks ];
           var21_x2y2 = Cv1v2[ node[l2]   + node[l1]   * num_blocks ];

           if(block_l1->getBlockLbPnt().getI() == 0 ) {
              var11_x12  = var12_x1y2 = 0;
           } else {
              var11_x12  = Cv1v1[ node_i[l1] + node[l2]   * num_blocks ];
              var12_x1y2 = Cv1v2[ node_i[l1] + node[l2]   * num_blocks ];
           }
             
           if(block_l2->getBlockLbPnt().getI() == 0 ) {
              var11_x21  = var21_x2y1 = 0;
           } else {
              var11_x21  = Cv1v1[ node[l1]   + node_i[l2] * num_blocks ];
              var21_x2y1 = Cv1v2[ node_i[l2] + node[l1]   * num_blocks ];
           }

           if(block_l1->getBlockLbPnt().getI() == 0 || 
              block_l2->getBlockLbPnt().getI() == 0 ) {
              var11_x1   = 0;
           } else {
              var11_x1   = Cv1v1[ node_i[l1] + node_i[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getJ() == 0 ) {
              var22_y12  = var21_x1y2 = 0;           
           } else {
              var22_y12  = Cv2v2[ node_j[l1] + node[l2]   * num_blocks ];
              var21_x1y2 = Cv1v2[ node[l2]   + node_j[l1] * num_blocks ];
           }
           
           if(block_l2->getBlockLbPnt().getJ() == 0 ) {
              var22_y21  = var12_x2y1 = 0;           
           } else {
              var22_y21  = Cv2v2[ node[l1]   + node_j[l2] * num_blocks ];
              var12_x2y1 = Cv1v2[ node[l1]   + node_j[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getJ() == 0 ||
              block_l2->getBlockLbPnt().getJ() == 0 ) {
              var22_y1   = 0;
           } else {
              var22_y1   = Cv2v2[ node_j[l1] + node_j[l2] * num_blocks ];
           }
           
           if(block_l1->getBlockLbPnt().getI() == 0 ||
              block_l2->getBlockLbPnt().getJ() == 0 ) {
              var12_x1y1 = 0.;  
           } else {
              var12_x1y1 = Cv1v2[ node_i[l1] + node_j[l2] * num_blocks ];
           }
           if(block_l1->getBlockLbPnt().getJ() == 0 ||
              block_l2->getBlockLbPnt().getI() == 0 ) {
              var21_x1y1 = 0.;
           }else {
              var21_x1y1 = Cv1v2[ node_i[l2] + node_j[l1] * num_blocks ];
           }
           
           C_vxvx = var11_x2   * x_x1[l1] * x_x1[l2]
                  + var11_x12  * x2_x[l1] * x_x1[l2]
                  + var11_x21  * x_x1[l1] * x2_x[l2]
                  + var11_x1   * x2_x[l1] * x2_x[l2];
           C_vyvy = var22_y2   * y_y1[l1] * y_y1[l2]
                  + var22_y12  * y2_y[l1] * y_y1[l2]
                  + var22_y21  * y_y1[l1] * y2_y[l2]
                  + var22_y1   * y2_y[l1] * y2_y[l2];
           C_vxvy = var12_x2y2 * x_x1[l1] * y_y1[l2]
                  + var12_x1y2 * x2_x[l1] * y_y1[l2] 
                  + var12_x2y1 * x_x1[l1] * y2_y[l2]
                  + var12_x1y1 * x2_x[l1] * y2_y[l2];
           C_vyvx = var21_x2y2 * y_y1[l1] * x_x1[l2]
                  + var21_x1y2 * y2_y[l1] * x_x1[l2]
                  + var21_x2y1 * y_y1[l1] * x2_x[l2]
                  + var21_x1y1 * y2_y[l1] * x2_x[l2];
           l1l2   = l1 + l2 * length;
           c_vx_vx[l1l2 ] = C_vxvx;
	   c_vy_vx[l1l2 ] = C_vyvx;
	   
	   if(       iflag_algo == 1 ) {
	      kernel[l1l2 ] = C_vyvy / vx[l1] / vx[l2];
	   } else if(iflag_algo == 2 ) {
              kernel[l1l2 ] = C_vxvx / vx[l1] / vx[l2]
                            + C_vyvy / vy[l1] / vy[l2]
		            - C_vxvy / vx[l1] / vy[l2]
                            - C_vyvx / vy[l1] / vx[l2];
	   } else {
	   } 
	   if(       iflag_algo == 1 ) {
	      f_new = C_vyvx / vx[l1];
	   } else if(iflag_algo == 2 ) {
   	      f_new = C_vyvx / vy[l1] - C_vxvx / vx[l1];	   
    	   } else {
	   } 

	   if(l1 == 0) {
	      c_eta_vx [l1l2] = 0;
	   } else {
	      c_eta_vx [l1l2] = c_eta_vx[l1 - 1 + l2 * length] 
		              + (f_new      + f_old    ) /2. * delta_x[l1];
	      f_old     = f_new;
	   }
       }
       if(iflag_algo == 2 ) {
          for(int l1 = 0; l1 < trajectory.size(); ++l1) {
	      l1l2   = l1 + l2 * length;    
              c_eta_vx[l1l2] *=	vy[l1]/ vx[l1];
          }
       }
   }
   // === ceta_dvx  ===
   double dvx_vx, dvy_vx;
   for(int l2 = 0; l2 < trajectory.size(); ++l2) {
       f_old = 0.;
       for(int l1 = 0; l1 < trajectory.size(); ++l1) {
	   int l1l2   = l1 +  l2      * length;
	   int l1l2m1 = l1 + (l2 - 1) * length;
	   int l1l2p1 = l1 + (l2 + 1) * length;

	   if(l2 == 0) {
	      if(fabs(delta_y[l2]) < 0.0000001) {
		 dvx_vx = 0.;
	         dvy_vx = 0.;	 
	      } else {
                 dvx_vx = (c_vx_vx[l1l2p1] - c_vx_vx[l1l2])/delta_y[l2];
	         dvy_vx = (c_vy_vx[l1l2p1] - c_vy_vx[l1l2])/delta_y[l2];	      
	      }
	   } else {
	      if(fabs(delta_y[l2]) < 0.0000001) {
		 dvx_vx = 0.;
	         dvy_vx = 0.;	 
	      } else {
                dvx_vx = (c_vx_vx[l1l2] - c_vx_vx[l1l2m1])/delta_y[l2];
	        dvy_vx = (c_vy_vx[l1l2] - c_vy_vx[l1l2m1])/delta_y[l2];
	      }
	   }

	   if(iflag_algo == 1 ) {
	      f_new = dvy_vx / vx[l1];
	   } else if(iflag_algo == 2 ) {
	      f_new = dvy_vx / vy[l1] - dvx_vx / vx[l1];
	   } else {
	   }

	   if(l1 == 0) {
	      c_eta_dvx [l1l2] = 0;
	   } else {
	      c_eta_dvx [l1l2] = c_eta_dvx[l1 - 1 + l2 * length] 
		              + (f_new      + f_old    ) /2. * delta_x[l1];
	      f_old     = f_new;
	   }
       }

       if(iflag_algo == 2 ) {
          for(int l1 = 0; l1 < trajectory.size(); ++l1) {
              int l1l2   = l1 +  l2      * length;
              c_eta_dvx[l1l2] *= 	vy[l1]/ vx[l1];
          }
       }
   }

   /*
   double dvx_vx, dvy_vx;
   for(int l2 = 0; l2 < trajectory.size(); ++l2) {
       f_old = 0.;
       for(int l1 = 0; l1 < trajectory.size(); ++l1) {
	   int l1l2   = l1 +  l2      * length;
	   int l1l2m1 = l1 + (l2 - 1) * length;
	   int l1l2p1 = l1 + (l2 + 1) * length;

	   if(l2 == 0) {
	      if(fabs(delta_y[l2]) < 0.0000001) {
		 dvx_vx = 0.;
	         dvy_vx = 0.;	 
	      } else {
                 dvx_vx = (c_vx_vx[l1l2p1] - c_vx_vx[l1l2])/delta_y[l2];
	         dvy_vx = (c_vy_vx[l1l2p1] - c_vy_vx[l1l2])/delta_y[l2];	      
	      }
	   } else {
	      if(fabs(delta_y[l2]) < 0.0000001) {
		 dvx_vx = 0.;
	         dvy_vx = 0.;	 
	      } else {
                dvx_vx = (c_vx_vx[l1l2] - c_vx_vx[l1l2m1])/delta_y[l2];
	        dvy_vx = (c_vy_vx[l1l2] - c_vy_vx[l1l2m1])/delta_y[l2];
	      }
	   }
	   
	   f_new = dvy_vx / vy[l1] - dvx_vx / vx[l1];
            
	   if(l1 == 0) {
	      c_eta_dvx [l1l2] = 0;
	   } else {
	      c_eta_dvx [l1l2] = c_eta_dvx[l1 - 1 + l2 * length] 
		              + (f_new      + f_old    ) /2. * delta_x[l1];
	      f_old     = f_new;
	   }
       }
       
       for(int l1 = 0; l1 < trajectory.size(); ++l1) {
           int l1l2   = l1 +  l2      * length;
           c_eta_dvx[l1l2] *= 	vy[l1]/ vx[l1];
       }
   }
   */
   // eta_eta	   
   int l1l2m1, l1m1l2, l1m1l2m1;
   for(int l2 = 0; l2 < length; ++l2) {
       l1l2 = l2 * length;
       c_eta_eta[l1l2] = 0.;
   }
   for(int l1 = 0; l1 < length; ++l1) {
       l1l2     = l1;
       c_eta_eta[l1l2] = 0.;
   }

   double tmp;
   for(int l2 = 1; l2 < length; ++l2) {
       for(int l1 = 1; l1 < length; ++l1) {
          l1l2     = l1     + l2 * length;   
          l1m1l2   = l1 - 1 + l2 * length;
          l1l2m1   = l1     + (l2 - 1) * length;   
          l1m1l2m1 = l1 - 1 + (l2 - 1) * length;
          tmp = ( kernel[l1l2  ] + kernel[l1m1l2] 
                + kernel[l1l2m1] + kernel[l1m1l2m1]
                )/4. * delta_x[l1] * delta_x[l2];
          c_eta_eta[l1l2] = c_eta_eta[l1m1l2] + c_eta_eta[l1l2m1] 
                          + tmp -  c_eta_eta[l1m1l2m1];
       }
   }

   if(iflag_algo == 2 ) {
      for(int l2 = 1; l2 < length; ++l2) {
          for(int l1 = 1; l1 < length; ++l1) {
              l1l2     = l1     + l2 * length;   
              c_eta_eta[l1l2] *= vy[l1] / vx[l1] * vy[l2] / vx[l2];
          }
      }
   }
   /*
   // ===== checking =====
   int l1l1, l2l2;
   for(int l2 = 0; l2 < length; ++l2) {
       for(int l1 = 0; l1 < length; ++l1) {
           l1l2     = l1 + l2 * length;
	   l1l1     = l1 + l1 * length;
	   l2l2     = l2 + l2 * length;
	   double stdv  = sqrt(c_eta_eta[l1l1] * c_eta_eta[l2l2] );
	   double stdev = sqrt(c_eta_eta[l1l1] *   c_vx_vx[l2l2] );
	   if(stdv > 0.00001 || stdev > 0.00001) {
               cout << l1 <<' ' << l2<<' ' << c_eta_eta[l1l2]/stdv << ' ' 
                                           << c_eta_vx[l1l2]/stdev << endl; 
	   }
       }
       cout<<endl;
   }
   exit(0);
   */
}

void Particle::calcTauAvg2Carte(double* c_eta_eta, double* c_eta_vx, 
		                double* c_eta_dvx, double* dvx_deta, double* c_vx_vx,
	                        double* delta_x, double* delta_y, double* vx){

   double tmp, tmp1, tmp2, tmp_i, tmp_im1;

   tau_Avg2Carte[0] = 0.;
   
   for(int i = 1; i < trajectory.size(); ++i) {
       int l   = i     +  i      * length;
       int lm1 = i - 1 + (i - 1) * length;
       
       if( i + 1 < trajectory.size() ) 
	   tmp_i = (dvx_deta[i + 1 ] - dvx_deta[i    ]) * 2./ (delta_y[i + 1 ] + delta_y[i     ]);
       if( i - 1 >= 0) 
	   tmp_im1 = (dvx_deta[i     ] - dvx_deta[i -1 ]) * 2./ (delta_y[i     ] + delta_y[i - 1 ]);
       if(i == 1                    ) tmp_im1 = tmp_i  ;
       if(i == trajectory.size() - 1) tmp_i   = tmp_im1;
	
       tmp_i = tmp_im1 = 0;     
       
       if(fabs(vx[i]) < 0.0000001) {
          tmp1 = 0.;
       } else {
	  tmp1 = 1. 
              - c_eta_dvx[ l ] / vx[i] 
	      - c_eta_eta[ l ] / 2. / vx[i] * tmp_i
	      + c_vx_vx[l]  / vx[i] / vx[i]
	      + 2. * c_eta_vx[l] / vx[i] / vx[i] * dvx_deta [i] 
	      + c_eta_eta[ l ] / vx[i] / vx[i] * dvx_deta [i] * dvx_deta [i]; 
          tmp1 /= vx[i]; 
       }
       
       if(fabs(vx[i - 1]) < 0.0000001) {
          tmp2 = 0.;
       } else {
	  tmp2 = 1.
	            - c_eta_dvx[ lm1 ] / vx[i-1] 
		    - c_eta_eta[ lm1 ] / 2. / vx[i-1] * tmp_im1
		    +   c_vx_vx[ lm1 ] / vx[i-1] / vx[i-1]
		    + 2. * c_eta_vx[ lm1 ] / vx[i-1] / vx[i-1] *  dvx_deta [i - 1]  
		    + c_eta_eta[ lm1 ] / vx[i-1] / vx[i-1] * dvx_deta [i-1] * dvx_deta [i-1]; 
	  tmp2  /= vx[i-1];
       }
       
       tmp  = (tmp1 + tmp2) * delta_x[i] /2.;	   
       tau_Avg2Carte[i] = tau_Avg2Carte[i - 1] + tmp;
   }
}

void Particle::calcTauVariCarte(double* vx, double* c_vx_vx, double* c_eta_eta,
                                double* dvx_deta, double* c_eta_vx, double* delta_x){
   double *f        = new double[length * length];
   
   int l1l2, l2l1;
   double v_avg2_l2, v_avg2_l1;
   for(int l2 = 0; l2 < length; l2++) {
      v_avg2_l2 = vx[l2] * vx[l2];
      for(int l1 = 0; l1 <  length; l1++) {
          v_avg2_l1 = vx[l1] * vx[l1];
	  
          l1l2 = l1 + l2 * length;
          l2l1 = l2 + l1 * length;
  
          f[l1l2] = (  c_vx_vx[   l1l2 ]
                     + c_eta_eta[ l1l2 ] * dvx_deta [l1] * dvx_deta [l2]
		     + c_eta_vx[  l1l2 ] * dvx_deta [l1]
		     + c_eta_vx[  l2l1 ] * dvx_deta [l2]
                    ) /v_avg2_l1 /v_avg2_l2;
      }
   }
   // --- now the 2-D integration (only need diagonal values) ---
   int l2l2;
   double tmp, var_old = 0.; 
   tau_VariCarte[0] = 0.0;
   for(int l2 = 1; l2 < length; l2++) {
      tmp = var_old;
      for(int l1 = 1; l1 < l2; l1++) {    // up row and right col
          l1l2 = l1 + l2 * length;
          tmp += (  f[l1l2 - 1] + f[l1l2 -     length]
                  + f[l1l2    ] + f[l1l2 - 1 - length]
                 ) / 4.0 
                  * delta_x[l1] * delta_x[l2] * 2.0;
      }
      l2l2 = l2 + l2 * length;
      tmp += (  f[l2l2] + f[l2l2 - 1] + f[l2l2 - 1] 
              + f[l2l2 - 1 - length] 
             )/4.0
              *  delta_x[l2] *  delta_x[l2];

      var_old = tmp;

      tau_VariCarte[l2] = tmp;
   }

   delete []f;
}

void Particle::PrintTauCarte(ostream & os) {
   if(try_Cartesian) {
      os.precision(9);
      os.setf(ios::fixed,ios::floatfield);
      os<<setw(12) <<tau_Avg0Carte[length-1] <<' ';
      os<<setw(12) <<tau_Avg2Carte[length-1] <<' ';
      os<<setw(12) <<tau_VariCarte[length-1] <<' ';
      os<<endl;
   } 	
}
// ===== the end of Cartesian system =====

// == 0th order ==
void Particle::calcTauAvg0(){
   double xzero = 1e-10;
   double vx, vx2, vx1, dvx, dx, x_len, dtx = 0.; 
   double vy, vy2, vy1, dvy, dy, y_len, dty = 0.;

   Block *block_i, *block_j;

   double t_sum = 0., dt = 0.;
   for(int i = 0; i < trajectory.size(); ++i) {
       block_i = domain_.point2block( trajectory[i] );
       x_len = trajectory[i].getX() 
             - block_i->getBlockLbPnt().getX();
       y_len = trajectory[i].getY() 
             - block_i->getBlockLbPnt().getY();
       dx    = block_i->getDx();
       dy    = block_i->getDy();
       vx1   = block_i->getVx1();
       vx2   = block_i->getVx2();
       vy1   = block_i->getVy1();
       vy2   = block_i->getVy2();
       
       dvx   = (vx2 - vx1 ) / dx;
       vx    = vx1 + dvx * x_len;
       dvy   = (vy2 - vy1 ) / dy;
       vy    =  vy1 + dvy * y_len;

       t_sum += dt;
       tau_Avg0[i] = t_sum;
        
       if(vx1 < 0 && vx < 0){ 
          if (fabs(dvx) < xzero) 
             dtx = (0. - x_len)/vx1;
          else 
             dtx = log(vx1/vx)/dvx;
       }
       if(vx2 > 0 && vx > 0){
         if (fabs(dvx) < xzero) 
             dtx = (dx - x_len)/vx2;
         else
             dtx = log(vx2/vx)/dvx;
       }
       if(vy1 < 0 && vy < 0){
         if (fabs(dvy) < xzero)
            dty = (0. - y_len)/vy1;
         else
            dty = log(vy1/vy)/dvy;
       }
       if(vy2 > 0 && vy > 0){
         if (fabs(dvy) < xzero) 
            dty = (dy - y_len)/vy2;
         else
            dty = log(vy2/vy)/dvy;
       }
       
       if(i + 1 < trajectory.size() ) {
          block_j = domain_.point2block(trajectory[i + 1]); 
          if(abs(  block_i->getBlockLbPnt().getI() 
                 - block_j->getBlockLbPnt().getI()
                ) == 0) { 
             dt = dty;
          }else {
             dt = dtx;
          }
       }
   }
}

// == 2nd order based on velocity components==
void Particle::calcTauAvg2(){
   double t_sum = 0., dt = 0.;
   double xzero = 1e-10;
   double vx, vx2, vx1, dvx, dx, x_len, x_len2, dtx = 0.; 
   double vy, vy2, vy1, dvy, dy, y_len, y_len2, dty = 0.;
   Block *block_i, *block_j;

   int num_intg = 50;
   double xtmp, ytmp;
   double * xx = new  double[num_intg];
   double * xw = new  double[num_intg];
   gauleg(-1.0, 1.0, xx, xw, num_intg);

   double x_bgn, y_bgn, x_end, y_end;
   double x1, x2, y1, y2;

   int num_blocks = domain_.getNumBlocks();
   double *Cv1v1  = domain_.getCv1v1();
   double *Cv2v2  = domain_.getCv2v2();
   double var11_x, var11_x1, var11_x2, var11_x12;
   double var22_y, var22_y1, var22_y2, var22_y12;
   int node_i1, node_j1, node_ij;

   for(int i = 0; i < trajectory.size(); ++i) {
       block_i = domain_.point2block( trajectory[i] );
  
       x_bgn   = trajectory[i].getX();
       y_bgn   = trajectory[i].getY(); 
       dx      = block_i->getDx();
       dy      = block_i->getDy();
       x1      = block_i->getBlockLbPnt().getX();
       y1      = block_i->getBlockLbPnt().getY();
       x2      = x1 + dx;
       y2      = y1 + dy;

       vx1   = block_i->getVx1();
       vx2   = block_i->getVx2();
       dvx   = (vx2 - vx1 ) / dx;
       vy1   = block_i->getVy1();
       vy2   = block_i->getVy2();
       dvy   = (vy2 - vy1 ) / dy;

       node_i1 = node_i[i];
       node_j1 = node_j[i];
       node_ij = node[i];
  
       //Pipat
	   int ii = block_i->getI();
	   int jj = block_i->getJ();
	   if((ii==0)||(jj==0)||(ii==domain_.getNx()-1)||(jj==domain_.getNy()-1))
	   {
		   var11_x1 = 0;
		   var11_x2 = 0;
           var11_x12= 0;
           var22_y1 = 0;
           var22_y2 = 0;
           var22_y12= 0;
	   }
	   else
	   {
		   var11_x1 = Cv1v1[ node_i1 + node_i1 * num_blocks ];
		   var11_x2 = Cv1v1[ node_ij + node_ij * num_blocks ];
           var11_x12= Cv1v1[ node_i1 + node_ij * num_blocks ];
           var22_y1 = Cv2v2[ node_j1 + node_j1 * num_blocks ];
           var22_y2 = Cv2v2[ node_ij + node_ij * num_blocks ];
           var22_y12= Cv2v2[ node_j1 + node_ij * num_blocks ];
	   }
	   //Pipat
	   /*Pipat
	   var11_x1 = Cv1v1[ node_i1 + node_i1 * num_blocks ];
       var11_x2 = Cv1v1[ node_ij + node_ij * num_blocks ];
       var11_x12= Cv1v1[ node_i1 + node_ij * num_blocks ];
       var22_y1 = Cv2v2[ node_j1 + node_j1 * num_blocks ];
       var22_y2 = Cv2v2[ node_ij + node_ij * num_blocks ];       
       var22_y12= Cv2v2[ node_j1 + node_ij * num_blocks ];
	   Pipat*/
       
       if(block_i->getBlockLbPnt().getI() == 0 )
          var11_x1 = var11_x12 = 0.;
       if(block_i->getBlockLbPnt().getJ() == 0 ) 
          var22_y1 = var22_y12 = 0.;
       
       t_sum += dt;
       tau_Avg2[i] = t_sum ;
       if(i + 1 < trajectory.size() ) {
          block_j = domain_.point2block(trajectory[i + 1]); 
          x_end   = trajectory[i + 1].getX();
          y_end   = trajectory[i + 1].getY();
          dty = 0;
          dtx = 0;
          for(int k = 0; k < num_intg; ++k) {
             xtmp  = (xx[k] * (x_end - x_bgn) + x2 + x_bgn)/2.;
             ytmp  = (xx[k] * (y_end - y_bgn) + y2 + y_bgn)/2.;
             x_len = xtmp - x1;
             y_len = ytmp - y1;
             x_len2= x2 - xtmp;
             y_len2= y2 - ytmp;
             vx    = vx1    + dvx    * x_len;
             vy    = vy1    + dvy    * y_len;

             var11_x =(var11_x2 * x_len  * x_len  
                     + var11_x1 * x_len2 * x_len2
                     + var11_x12* x_len  * x_len2 * 2.)/dx/dx;
             
             var22_y =(var22_y2 * y_len  * y_len 
                     + var22_y1 * y_len2 * y_len2 
                     + var22_y12* y_len  * y_len2 * 2.)/dy/dy;

             dtx  += xw[k] * 1./vx *(1. + var11_x/vx/vx); 
             dty  += xw[k] * 1./vy *(1. + var22_y/vy/vy);
          }
          if(abs(  block_i->getBlockLbPnt().getI() 
                 - block_j->getBlockLbPnt().getI()
                ) == 0
            ) { 
             dt = dty * (y_end - y_bgn)/2.;
          }else {
             dt = dtx * (x_end - x_bgn)/2.;
          }
       }
   }
   delete []xx;
   delete []xw;
}

void Particle::calcTauAvg2Trans(double *Eta_Eta, double *Vxi_Eta) {

   double* etadu = new double [trajectory.size()];
   calcEtaDU(etadu);

   double t_sum = 0.;
   double dt    = 0.;
   double tmp1, tmp2;

   // 0th-order
   tau_Avg0Trans[0] = t_sum;
   for(int i = 1; i < trajectory.size(); ++i) {
       tmp1   = 1./ v_avg[i  ];
       tmp2   = 1./ v_avg[i-1];
       dt     = 0.5 * (tmp1 + tmp2) * ds[i];
       t_sum += dt;
       tau_Avg0Trans[i] = t_sum;
   }

   // 2th-order
   t_sum = 0.;
   tau_Avg2Trans[0] = t_sum;
   for(int i = 1; i < trajectory.size(); ++i) {
       tmp1   = 1./ v_avg[i  ] 
	      * (1. 
	            + v_var[i]/ v_avg[i]/v_avg[i]
		 //- etadu[i] //The correct derivation shouldn't have this term (Pipat)
		    + Eta_Eta[i + i * length] / v_avg[i]/v_avg[i] * dv_deta[i] * dv_deta[i]
		    + 2. * Vxi_Eta[i *(1 + length)] * dv_deta[i] / v_avg[i] / v_avg[i]
                );
       tmp2   = 1./ v_avg[i-1] 
	      * (1. 
		    + v_var[i-1] / v_avg[i-1]/v_avg[i-1]
		 //- etadu[i-1] //The correct derivation shouldn't have this term (Pipat)
		    + Eta_Eta[(i - 1) * (1+ length)] / v_avg[i-1]/v_avg[i-1] * dv_deta[i-1] * dv_deta[i-1]
		    + 2. * Vxi_Eta[(i - 1) *(1+length)] * dv_deta[i-1] / v_avg[i-1] / v_avg[i-1]
                );                       
       dt     = 0.5 * (tmp1 + tmp2) * ds[i];
       t_sum += dt;
       tau_Avg2Trans[i] = t_sum;
   }

   delete []etadu;
}

// == TauVar ==
void Particle::calcTauVar(){
   int l1l2, l2l1;
   double v_avg2_l1, v_avg2_l2;

   double *f        = new double[length * length];
   double *Vxi_Vxi  = new double[length * length];
   double *Vxi_Eta  = new double[length * length];
   double *Eta_Eta  = new double[length * length];
   
   CoordTrans2(Vxi_Vxi, Vxi_Eta, Eta_Eta); 
   calcTauAvg2Trans(Eta_Eta, Vxi_Eta);
   
   for(int l2 = 0; l2 < length; l2++) {
      v_avg2_l2 = v_avg[l2]*v_avg[l2];
      for(int l1 = 0; l1 <  length; l1++) {
          v_avg2_l1 = v_avg[l1] * v_avg[l1];
          l1l2 = l1 + l2 * length;
          l2l1 = l2 + l1 * length;
          f[l1l2] = ( Vxi_Vxi[l1l2]
                    + Vxi_Eta[l1l2] * dv_deta[l2]
                    + Vxi_Eta[l2l1] * dv_deta[l1]
                    + Eta_Eta[l1l2] * dv_deta[l1] * dv_deta[l2]
                    ) /v_avg2_l1 /v_avg2_l2;
	  //cout<<l1<<" "<<l2<<" "<<f[l1l2]<<" "<<Vxi_Vxi[l1l2]<<" "<<Eta_Eta[l1l2]<<" "<<Vxi_Eta[l1l2]<<" "<<Vxi_Eta[l2l1]<<" "<<dv_deta[l2]<<" "<<dv_deta[l1]<<endl; 
      }
   }
   delete [] Vxi_Vxi;
   delete [] Vxi_Eta;
   delete [] Eta_Eta;
   
   // --- now the 2-D integration (only need diagonal values) ---
   int l2l2;
   double tmp, var_old = 0.; 
   tau_VariTrans[0] = 0.0;
   for(int l2 = 1; l2 < length; l2++) {
      tmp = var_old;
      for(int l1 = 1; l1 < l2; l1++) {    // up row and right col
          l1l2 = l1 + l2 * length;
          tmp += (  f[l1l2 - 1] + f[l1l2 - length]
                  + f[l1l2] + f[l1l2 - 1 - length]
                 ) / 4.0 
                  * ds[l1] * ds[l2] * 2.0;
      }
      l2l2 = l2 + l2 * length;    
      tmp += (  f[l2l2] + f[l2l2 - 1] + f[l2l2 - 1] 
              + f[l2l2 - 1 - length] 
             )/4.0
              * ds[l2] * ds[l2];

      var_old = tmp;

      tau_VariTrans[l2] = tmp;
   }
   delete []f;
}
// ===== the end of transformed system =====

// ===== OUPPUT =====
void Particle::PrintParticle(ostream & os){
   os.precision(9);
   os.setf(ios::fixed,ios::floatfield);
   os<<setw(12)<<particle_point.getX()<<' ';
   os<<setw(12)<<particle_point.getY()<<' ';
   os<<setw(12)<<particle_point.getI()<<' ';
   os<<setw(12)<<particle_point.getJ()<<' ';
   os<<setw(12)<<particle_time_       <<' ';
   os<<endl;
}

void Particle::PrintTrajectory(ostream & os){
   os.precision(9);
   os.setf(ios::fixed,ios::floatfield);
   os<<setw(12) <<tau_Avg0Trans[length-1] <<' ';
   os<<setw(12) <<tau_Avg2Trans[length-1] <<' ';
   os<<setw(12) <<tau_VariTrans[length-1] <<' ';
   os<<setw(12) <<tau_Avg0[length-1]      <<' ';
   os<<setw(12) <<tau_Avg2[length-1]      <<' ';
   os<<setw(12) <<tau_Vari[length-1]      <<' ';
   os<<setw(12) <<tau_AvgPdf[length-1]    <<' ';
   os<<endl;
}

void Particle::PrintPathline(ostream & os){
   os.precision(9);
   double dist, tmp, xp,xpo, yp, ypo;
    
   os.setf(ios::fixed,ios::floatfield);
   for(int i = 0; i < trajectory.size();++i) {
      xp = trajectory[i].getX();
      yp = trajectory[i].getY();
      if(i == 0) {
	 dist = 0;
      } else {
	 xpo = trajectory[i-1].getX();
         ypo = trajectory[i-1].getY();
	 tmp =  sqrt(  (xp - xpo ) * ( xp - xpo ) + (yp - ypo ) * ( yp - ypo ) );
	 dist += tmp;
      }
      os<<setw(12) << xp   << ' ';
      os<<setw(12) << yp   << ' ';
      os<<setw(12) << dist << ' ';
      
      os<<setw(12) << tau_Avg0Trans[i]     <<' ';
      os<<setw(12) << tau_Avg2Trans[i]     <<' ';
      os<<setw(12) << tau_VariTrans[i]     <<' ';
      os<<setw(12) << tau_Avg0Carte[i]     <<' ';
      os<<setw(12) << tau_Avg2Carte[i]     <<' ';
      os<<setw(12) << tau_VariCarte[i]     <<' ';
      os<<setw(12) << tau_Avg0[i]          <<' ';
      os<<setw(12) << tau_Avg2[i]          <<' ';
      os<<setw(12) << tau_Vari[i]          <<' ';
      

    //os<<setw(12) << time_flight[i]       <<' ';
      os<<endl;
   }
   os<<endl;
}

void Particle::PrintVelocity(ostream & os){
   os.precision(9);
   os.setf(ios::fixed,ios::floatfield);
   os<<" Velo_Avg     Velo_Var   " << endl;
   for(int i = 0; i < trajectory.size();++i) {
      os<<setw(12)<<v_avg[i] << ' ';
      os<<setw(12)<<v_var[i] << endl;
   }
   os<<endl;
}

void Particle::PrintTrajectory(){
   for(int i = 0; i < trajectory.size();++i) {
      trajectory[i].PrintXY(cout);
   }
}

void Particle::calcTauAvgPdf(){
   int num_intg = 50;
   double * xx = new double[num_intg];
   double * xw = new double[num_intg];
   gauleg(-1.0, 1.0, xx, xw, num_intg);

   double t_sum = 0.;
   double dt    = 0.;
   double tmp1, tmp2, xtmp;
   double vtmp1, vtmp2;
   tau_AvgPdf[0] = t_sum;
   for(int i = 1; i < trajectory.size(); ++i) {
       if(v_var[i-1] < 0.0001) {
          tmp1 = 1./ v_avg[i-1];
       } else {
          vtmp1 = 0.0001 * v_avg[i-1];     
          vtmp2 = v_avg[i-1] + 500. * sqrt(v_var[i-1]); 
          tmp1 = 0.;
          for(int kt = 0; kt < num_intg; ++kt) {
              xtmp  = (xx[kt] * (vtmp2 - vtmp1) + vtmp2 + vtmp1)/2.;
              tmp1 += 1./xtmp * pdfLognormal(xtmp,v_avg[i-1],v_var[i-1]) * xw[kt];  
          }
          tmp1 *= (vtmp2 - vtmp1)/2.;
       }
       if(v_var[i] < 0.0001) {
          tmp2 = 1./ v_avg[i];
       } else {
          vtmp1 = 0.0001 * v_avg[i];     
          vtmp2 = v_avg[i] + 500. * sqrt(v_var[i]); 
          tmp2 = 0.;
          for(int kt = 0; kt < num_intg; ++kt) {
              xtmp  = (xx[kt] * (vtmp2 - vtmp1) + vtmp2 + vtmp1)/2.;
              tmp2 += 1./xtmp * pdfLognormal(xtmp,v_avg[i],v_var[i]) * xw[kt];  
          }
          tmp2 *= (vtmp2 - vtmp1)/2.;
       }
       dt     = 0.5 * (tmp1 + tmp2) * ds[i];
       t_sum += dt;
       tau_AvgPdf[i] = t_sum;
   }
   delete []xw;
   delete []xx;
}

