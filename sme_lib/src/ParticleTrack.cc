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

#include "ParticleTrack.h"

ParticleTrack::ParticleTrack(Domain* domain, int nstrl, Point* pt, 
                      bool print_option)
:  domainPtr(domain), nstrl_(nstrl), print_status(print_option)
{  
   debug = false;

   unit = domainPtr->getUnit();
   
   particles = new Particle* [nstrl_];
   for (int i = 0; i < nstrl_; ++i) {
       particles[i] =  new Particle( pt[i], 0.0, *domainPtr);
   }
   initialization();
}

// the following one was used!
ParticleTrack::ParticleTrack(Domain* domain, bool debug_, bool print_option)
:  domainPtr(domain), print_status(print_option), debug(debug_)
{
   // debug = true;	
   cout << "ParticleTrack::ParticleTrack()" << endl;	
   unit = domainPtr->getUnit();
   nstrl_ = 0;

   // 
   for(int i = 0; i < domainPtr->getNumWells(); ++i) {
       if( domainPtr->getWellArc(i)->isProdWell() ) {	   
           prodWellArcPtr = domainPtr->getWellArc(i); 
       }
       if( domainPtr->getWellArc(i)->isInjcWell() ) {
           nstrl_ = nstrl_ + domainPtr->getWellArc(i)->length();
       }
   }
   for(int i = 0; i < 6; i++) {
       if( domainPtr->getContr()->isInjcBound(i) ) {
           nstrl_ = nstrl_ + domainPtr->getContr()->getBoundLength(i);
       }
   }
   particles = new Particle* [nstrl_];
   int nstrl_old = 0;
   for(int i = 0; i < domainPtr->getNumWells(); ++i) {
       if( domainPtr->getWellArc(i)->isInjcWell() ) {
           int nstrl_well = domainPtr->getWellArc(i)->length();
           for(int j = 0; j < nstrl_well; ++j) {
	       //cout << "i,j:" << i <<' ' << j <<' '<< j + nstrl_old << endl;
               particles[j + nstrl_old] = new Particle(
                      (*domainPtr->getWellArc(i))[j], 0.0, *domainPtr
                                                      );
           }
           nstrl_old += nstrl_well;
       }
   }
   for (int i = 0; i < 6; i++) {
       if ( domainPtr->getContr()->isInjcBound(i) ) {
           int nstrl_boundary = domainPtr->getContr()->getBoundLength(i);
           for (int j = 0; j < nstrl_boundary; j++) {
               particles[j + nstrl_old] = new Particle( 
                       domainPtr->getContr()->getBoundPts(i)[j], 0.0, *domainPtr
                                                       );  
           }
           nstrl_old += nstrl_boundary; 
       }
   } 
   //exit(0);
   initialization();
}

ParticleTrack::~ParticleTrack(){
   for(int i = 0; i < nstrl_; ++i)
       delete particles [i];
   delete [] particles;
   
   delete [] trvl_avg;
   delete [] trvl_var;
   delete [] trvl_rho;
   delete [] trvl_avg_log;
   delete [] trvl_var_log;
   delete [] trvl_rho_log;
   delete [] vbar;
   delete [] Cqq;
}

void::ParticleTrack::initialization() {
    trvl_avg     = new double [ nstrl_ ];
    trvl_var     = new double [ nstrl_ ];
    trvl_rho     = new double [ nstrl_ * nstrl_ ];
    trvl_avg_log = new double [ nstrl_ ];
    trvl_var_log = new double [ nstrl_ ];
    trvl_rho_log = new double [ nstrl_ * nstrl_];
    vbar         = new double [ nstrl_ ];
    Cqq          = new double [ nstrl_ * nstrl_];

    Cv1v1   = domainPtr->getCv1v1();
    Cv2v2   = domainPtr->getCv2v2();
    Cv1v2   = domainPtr->getCv1v2();
}

void ParticleTrack::goTrack() {
   if(debug) cout << "ParticleTrack Starts " << endl;
   int num_strl_prod = 0;

   int * iwp = new int [domainPtr->getNumWells()];
   for (int iw = 0; iw < domainPtr->getNumWells(); ++iw)
        iwp[iw] = 0;
   
   for(int i = 0; i < nstrl_; ++i) {
       if(debug) cout<<"Particle "<< i <<" is moving..."<<endl;
       
       while( !particles[i]->getBlock()->isProdWell() && !particles[i]->getBlock()->isProdBound() ) {
          particles[i]->move();
       }
       //cout << i <<' ' << particles[i]->getBlock()->isProdWell() << endl;
       
       if(particles[i]->getBlock()->isProdWell()) {

			for(int iw = 0; iw < domainPtr->getNumWells(); ++iw) {
				if(domainPtr->getWellArc(iw)->isProdWell() ) {
					if( particles[i]->getBlock()->getI() ==  domainPtr->getWellArc(iw)->getI() &&
						particles[i]->getBlock()->getJ() ==  domainPtr->getWellArc(iw)->getJ()){
		     /*	 
		     cout << i<< ' '
			  << domainPtr->getWellArc(iw)->getI() << ' '
		          << domainPtr->getWellArc(iw)->getJ() << ' '	  
                          << iw << ' ' 
			  << domainPtr->getWellArc(iw)->getId() << endl;
			 */
						++iwp[iw];
						int id = iwp[iw] - 1;
						domainPtr->getWellArc(iw)->setPoint(id, particles[i]->getBlock()->getI(),
		  					particles[i]->getBlock()->getJ(),
				  			particles[i]->getTrajectory().getX(),
					   		particles[i]->getTrajectory().getY());
					}
				}
			}
	   }
       
       particles[i]->CoordTrans();
       particles[i]->calcTauAvg();
       particles[i]->calcTauVar();
   } 

   // iwp[iw]: local variabla= the number of particles at well iw.
   for(int iw = 0; iw < domainPtr->getNumWells(); ++iw) {
       if(domainPtr->getWellArc(iw)->isProdWell() ) {
			domainPtr->getWellArc(iw)->setLength(iwp[iw]);
			domainPtr->getWellArc(iw)->setPointOrder( domainPtr->getNx(),
					            domainPtr->getNy());
	  //cout << iwp[iw] << endl;
	  //domainPtr->getWellArc(iw)->DisplayArc();
       }
   }
   delete [] iwp;
   // exit(0);

   // print out
   if(print_status) {

      ofstream fileOu ("pathline.out", ios::out);
      /*
      cout << "in debugging.." << endl;
      for(int i = 0; i < nstrl_; ++i){
          particles[i]->PrintPathline(fileOu);
      }
      fileOu.close();
      exit(0);*/

      ofstream fileOu1("Vbar_AlongPath.out", ios::out);
      ofstream fileOu2("TravelTimes.out", ios::out);
      ofstream fileOu3("TravelTimes_Carte.out", ios::out);
      if ( (!fileOu) || (!fileOu1) ){
         cerr << "\n File does not exist !\n";
         exit(EXIT_FAILURE);
      }
      
      for(int i = 0; i < nstrl_; ++i){
          particles[i]->PrintPathline(fileOu);
	  particles[i]->PrintVelocity(fileOu1);
	  particles[i]->PrintTrajectory(fileOu2);
	  particles[i]->PrintTauCarte(fileOu3);
      }
      fileOu.close();
      fileOu1.close();
      fileOu2.close();
      fileOu3.close();
   }
}

/* Function: goTrackBackward
 * Author: Pipat Likanapaisal
 * --------------------------
 * The function utilizes Particle::moveBackward. Then, the
 * functions call Particle::reorderData to arrange trajectory
 * and time_flight in forward tracking fashion so that the
 * functions for statistical calculation can be used.
 */
void ParticleTrack::goTrackBackward() 
{
	if(debug) cout << "ParticleBackTrack Starts " << endl;
	for(int ii=0;ii<nstrl_;ii++)
	{
		if(debug) cout<<"Particle "<< ii <<" is moving..."<<endl;
		while(!particles[ii]->getBlock()->isInjcWell())
		{
			particles[ii]->moveBackward();
		}
		particles[ii]->reorderData();
		particles[ii]->CoordTrans();
		particles[ii]->calcTauAvg();
		particles[ii]->calcTauVar();
	}
	if(print_status) {
		ofstream fileOu ("pathline.out", ios::out);
		ofstream fileOu1("Vbar_AlongPath.out", ios::out);
		ofstream fileOu2("TravelTimes.out", ios::out);
		ofstream fileOu3("TravelTimes_Carte.out", ios::out);
		if ( (!fileOu) || (!fileOu1) )
		{
			cerr << "\n File does not exist !\n";
			exit(EXIT_FAILURE);
		}
		for(int i = 0; i < nstrl_; ++i)
		{
			particles[i]->PrintPathline(fileOu);
			particles[i]->PrintVelocity(fileOu1);
			particles[i]->PrintTrajectory(fileOu2);
			particles[i]->PrintTauCarte(fileOu3);
		}
		fileOu.close();
		fileOu1.close();
		fileOu2.close();
		fileOu3.close();
	}
}

void ParticleTrack::calcTrvlTimeMoments() {
   cout <<"ParticleTrack::calcTrvlTimeMoments() " << endl;
   for(int iw = 0; iw < domainPtr->getNumWells(); ++iw) {
       cout << "iw = " << iw << endl;
       if(domainPtr->getWellArc(iw)->isProdWell() ) {

			prodWellArcPtr = domainPtr->getWellArc(iw); 
			prodBlkPtr     = domainPtr->getWellArc(iw)->getBlock();

			int i_prod  = prodBlkPtr->getI();
			int j_prod  = prodBlkPtr->getJ();
			int k_prod  = prodBlkPtr->getK();

			double poro = domainPtr->getPoro(i_prod, j_prod, k_prod);
			double cell_x_len = prodBlkPtr->getDx();
			double cell_y_len = prodBlkPtr->getDy();
			double cell_z_len = prodBlkPtr->getDz();

			double vbar_total;
			vbar_total = (   fabs( prodBlkPtr->getVx2() - prodBlkPtr->getVx1() ) * cell_y_len   
							+ fabs( prodBlkPtr->getVy2() - prodBlkPtr->getVy1() ) * cell_x_len
						) * cell_z_len * poro;

			if(unit == FIELD ) {
				cout << "Total Production = " <<  vbar_total / stbPday_m3Psec << " stb/day " << endl;
			}
			else {
				cout << "Total Production = " <<  vbar_total << endl;
			}

			double * weight = prodWellArcPtr->getWeight();
			double   radius = prodWellArcPtr->getRadius();
			double * theta  = prodWellArcPtr->getTheta();
			double * v_arc  = prodWellArcPtr->getTotalV();
	  
			double vbar_total_int = 0.;
			for(int i = 0; i < prodWellArcPtr->length(); ++i) {
				vbar[i] = v_arc[i] * theta[i] * weight[i] * poro;
				//cout << v_arc[i] << ' ' << weight[i]<<' ' << poro << endl;
				vbar_total_int += vbar[i];
			}
			if(unit == FIELD ) {
				cout << "Total Production = " <<  vbar_total_int / stbPday_m3Psec << " stb/day " << endl;
			}
			else {
				cout << "Total Production = " << vbar_total_int << endl;
			}

			double sum = 0.;
			for(int i = 0; i < prodWellArcPtr->length(); ++i) {
				vbar[i] = vbar[i] * vbar_total / vbar_total_int;
				sum += vbar[i];
			}

			if(unit == FIELD )
				cout << "Total Injection = " <<   sum / stbPday_m3Psec << " by Field Unit" << endl;
			else
				cout << "Total Injection = " <<  vbar_total << endl;

			//exit(0);
			for(int i = 0; i < nstrl_ ; ++i) {
				trvl_avg[i] = particles[i]->getLastTravelTime();
				trvl_var[i] = particles[i]->getLastTravelTVar();
			}
			transform1(nstrl_, trvl_avg, trvl_var, trvl_avg_log, trvl_var_log);
	  
			calcTrvlTimeCorrelation();

			transform2(nstrl_, trvl_avg, trvl_var_log, trvl_rho, trvl_rho_log);

			calcCqqAtWell(radius, weight);
			for(int p2 = 0; p2 < nstrl_; ++p2){
				for(int p1 = 0; p1 < nstrl_; ++p1) {
					Cqq[p1 + p2 * nstrl_ ] *= poro * poro;
				}
			}
       }
   }
	
   //exit(0);
}

void ParticleTrack::calcTrvlTimeMomentsNoWell() {

  for(int i = 0; i < nstrl_ ; ++i) {
      trvl_avg[i] = particles[i]->getLastTravelTime();
      trvl_var[i] = particles[i]->getLastTravelTVar();
  }

  
  transform1(nstrl_, trvl_avg, trvl_var, trvl_avg_log, trvl_var_log);
  
  calcTrvlTimeCorrelation();

  transform2(nstrl_, trvl_avg, trvl_var_log, trvl_rho, trvl_rho_log);
}

void ParticleTrack::calcTrvlTimeCorrelation() {
   int p1p2, l1l2, l2l1, l1m1l2, l1l2m1, l1m1l2m1;

   int len2, len1;
   double v_avg1, v_avg2, v_avg11, v_avg22;
   double ds1, ds2;
   double dvdeta2, dvdeta1;
   double sum;
   
   int num_max = 0;
   for(int p = 0; p < nstrl_; ++p){
       if( particles[p]->getLength() > num_max ) 
           num_max = particles[p]->getLength();
   }
   double *f         = new double [num_max * num_max];
   double *Vxi_Vxi   = new double [num_max * num_max];
   double *Vxi_Veta  = new double [num_max * num_max];
   double *Eta_Eta   = new double [num_max * num_max];
   double *Vxi_Eta   = new double [num_max * num_max];
   double *Vxi_Eta2  = new double [num_max * num_max];
   double *Veta_Veta = new double [num_max * num_max];
       
   for(int p2 = 0; p2 < nstrl_; ++p2){
       len2 = particles[p2]->getLength();
       for(int p1 = 0; p1 < nstrl_; ++p1) {
          len1 = particles[p1]->getLength(); 
          for(int l = 0; l < num_max * num_max; ++l ) {
                     f[l] = 0.;
               Vxi_Vxi[l] = 0.;
              Vxi_Veta[l] = 0.;
               Eta_Eta[l] = 0.;
               Vxi_Eta[l] = 0.;
              Vxi_Eta2[l] = 0.;
             Veta_Veta[l] = 0.;
          }
          calcMoments(particles[p2], particles[p1], Vxi_Veta, Vxi_Eta2);
          calcMoments(particles[p1], particles[p2], Vxi_Vxi,  Vxi_Eta,
                      Eta_Eta, Vxi_Veta, Veta_Veta);

	  //Due to the fact that dv_deta is cell-center value, f[] will be
	  //calculated at the time of integration.
	  /*
          for(int l2 = 0; l2 < len2; l2++) {
              v_avg2 = particles[p2]->getVAvg(l2);          
              v_avg22= v_avg2 * v_avg2;
                 ds2 = particles[p2]->getDs(l2);   
             dvdeta2 = particles[p2]->getDvDeta(l2);
              for(int l1 = 0; l1 < len1; l1++) {
                 v_avg1 = particles[p1]->getVAvg(l1);
                 v_avg11= v_avg1 * v_avg1;
                    ds1 = particles[p1]->getDs(l1);
                dvdeta1 = particles[p1]->getDvDeta(l1);
                 l1l2 = l1 + l2 * len1;
                 l2l1 = l2 + l1 * len2;
                 f[l1l2] = ( Vxi_Vxi[l1l2]
                           + Vxi_Eta[l1l2] * dvdeta2
                           + Vxi_Eta2[l2l1]* dvdeta1
                           + Eta_Eta[l1l2] * dvdeta1 * dvdeta2
                           ) /v_avg22 / v_avg11;
              }
          }
	  */

	  int l1m[4] = {0,-1,0,-1};
	  int l2m[4] = {0,0,-1,-1};
          sum = 0.;
          for(int l2 = 1; l2 < len2; l2++) {
	      dvdeta2 = particles[p2]->getDvDeta(l2);
              for(int l1 = 1; l1 < len1; l1++) {
		  /*
		  l1l2     = l1     +  l2      * len1;
		  l1m1l2   = l1 - 1 +  l2      * len1;
		  l1l2m1   = l1     + (l2 - 1) * len1;   
		  l1m1l2m1 = l1 - 1 + (l2 - 1) * len1;
		  sum +=(  f[l1l2  ] + f[l1m1l2]  
		  	 + f[l1l2m1] + f[l1m1l2m1] 
			)/4. * particles[p1]->getDs(l1)
		             * particles[p2]->getDs(l2);
		  */

		  double f_sum = 0;
		  dvdeta1 = particles[p1]->getDvDeta(l1);
		  for(int pp=0;pp<4;pp++){
		    l1l2 = (l1+l1m[pp])+(l2+l2m[pp])*len1;
		    l2l1 = (l2+l2m[pp])+(l1+l1m[pp])*len2;
		    f_sum += ( Vxi_Vxi[l1l2]
			       + Vxi_Eta[l1l2] * dvdeta2
			       + Vxi_Eta2[l2l1]* dvdeta1
			       + Eta_Eta[l1l2] * dvdeta1 * dvdeta2
			     ) /pow(particles[p2]->getVAvg(l2+l2m[pp]),2.) 
		               /pow(particles[p1]->getVAvg(l1+l1m[pp]),2.);
		  }
		  sum += f_sum/4.*particles[p1]->getDs(l1)*particles[p2]->getDs(l2);
              }
          }
          p1p2 = p1 + p2 * nstrl_ ;
          trvl_rho[p1p2] = sum;
       }
   }
   /*
   // Symmetric Checking!!
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = 0; p1 < nstrl_; ++p1) {
           double error = fabs( trvl_rho[p1 + p2 * nstrl_] - trvl_rho[p2 + p1 * nstrl_] );
           double value1 = fabs( trvl_rho[p1 + p2 * nstrl_]);
           double value2 = fabs( trvl_rho[p2 + p1 * nstrl_]);
           double value = value1;
           if( error > 1.0e-8) {
                 if(value2 < value) value = value2;
               double error_rel = error/value;
               if(error_rel > 1.0e-8 ) {
                  cout << p1 << ' '
                    << p2 << ' '
                    << trvl_rho[p1 + p2 * nstrl_] << ' '
                    << trvl_rho[p2 + p1 * nstrl_] << ' '
                    << trvl_rho[p1 + p2 * nstrl_] - trvl_rho[p2 + p1 * nstrl_]
                    << endl;
               }
           }
       }
   }
   //exit(0);
   */
   delete[] Vxi_Eta;
   delete[] Vxi_Eta2;
   delete[] Eta_Eta;
   delete[] Vxi_Veta;
   delete[] Vxi_Vxi;
   delete[] Veta_Veta; 
   delete[] f;

   // Symmetry
   /*
   int p2p1;
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = p2 + 1; p1 < nstrl_; ++p1){
          p1p2 = p1 + p2 * nstrl_ ;
          p2p1 = p2 + p1 * nstrl_ ;
          trvl_rho[p1p2] = trvl_rho[p2p1];
       }
   }
   */
   ofstream fileOu("TravelTimeCorr.out", ios::out);
   fileOu << setiosflags(ios::fixed | ios::showpoint) << setprecision(3);
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = 0; p1 < nstrl_; ++p1){        
           p1p2 = p1 + p2 * nstrl_ ;    
           fileOu<< p1 <<' ' 
		 << p2 <<' '
		 << trvl_rho[p1p2]<< endl;
       }
       fileOu<<endl;
   }
   fileOu.close();
   //cout <<"stop at ParticleTracking.out " << endl;
   //exit(0);
}

void ParticleTrack::calcCqqAtWell(double radius, double* weight) {
   Block *block_l1, *block_l2;
   double C_vxvx, C_vyvy, C_vxvy, C_vyvx;
   int node_l1, node_l2;
   int node_i_l1, node_i_l2;
   int node_j_l1, node_j_l2;
   double var11_x , var11_x1,   var11_x2,   var11_x12,  var11_x21;
   double var22_y , var22_y1,   var22_y2,   var22_y12,  var22_y21;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double var21_xy, var21_x1y1, var21_x1y2, var21_x2y1, var21_x2y2;
   int num_blocks  = domainPtr->getNumBlocks();
   int l1l2;
   double  costh_l2, sinth_l2, costh_l1, sinth_l1;
   double  x_x1_l2, x2_x_l2, y_y1_l2, y2_y_l2;
   double  x_x1_l1, x2_x_l1, y_y1_l1, y2_y_l1;
   
   /*
   for(int p2 = 0; p2 < nstrl_; ++p2){
        cout << particles[p2]->getTrajectory(0).getX() << ' '
             << particles[p2]->getTrajectory(0).getY() << endl;
   }
   exit(0);
   */
   for(int p2 = 0; p2 < nstrl_; ++p2){
        block_l2 = domainPtr->point2block( particles[p2]->getTrajectory(0) );
         node_l2 = particles[p2]->getNode(0);
       node_i_l2 = particles[p2]->getNode_i(0);
       node_j_l2 = particles[p2]->getNode_j(0);
        costh_l2 = particles[p2]->getCosth(0);
        sinth_l2 = particles[p2]->getSinth(0);
         x_x1_l2 = particles[p2]->getX_x1(0);
         x2_x_l2 = particles[p2]->getX2_x(0);
         y_y1_l2 = particles[p2]->getY_y1(0);
         y2_y_l2 = particles[p2]->getY2_y(0);

       for(int p1 = 0; p1 < nstrl_; ++p1) {
            block_l1 = domainPtr->point2block( particles[p1]->getTrajectory(0) );
             node_l1 = particles[p1]->getNode(0);
           node_i_l1 = particles[p1]->getNode_i(0);
           node_j_l1 = particles[p1]->getNode_j(0);
            costh_l1 = particles[p1]->getCosth(0);
            sinth_l1 = particles[p1]->getSinth(0);
             x_x1_l1 = particles[p1]->getX_x1(0);
             x2_x_l1 = particles[p1]->getX2_x(0);
             y_y1_l1 = particles[p1]->getY_y1(0);
             y2_y_l1 = particles[p1]->getY2_y(0);
              
          var11_x2   = Cv1v1[ node_l1   + node_l2   * num_blocks ];
          var22_y2   = Cv2v2[ node_l1   + node_l2   * num_blocks ];
          var12_x2y2 = Cv1v2[ node_l1   + node_l2   * num_blocks ];
          var21_x2y2 = Cv1v2[ node_l2   + node_l1   * num_blocks ];
          if(block_l1->getBlockLbPnt().getI() == 0 ) {
             var11_x12  = var12_x1y2 = 0;
          } else { 
             var11_x12  = Cv1v1[ node_i_l1 + node_l2  * num_blocks ];
             var12_x1y2 = Cv1v2[ node_i_l1 + node_l2  * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x21  = var21_x2y1 = 0;
          } else {
             var11_x21  = Cv1v1[ node_l1   + node_i_l2 * num_blocks ];
             var21_x2y1 = Cv1v2[ node_i_l2 + node_l1   * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ) {
             var22_y12  = var21_x1y2 = 0;          
          } else {
             var22_y12  = Cv2v2[ node_j_l1 + node_l2   * num_blocks ];
             var21_x1y2 = Cv1v2[ node_l2   + node_j_l1 * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y21  = var12_x2y1 = 0;          
          } else {
             var22_y21  = Cv2v2[ node_l1   + node_j_l2 * num_blocks ];
             var12_x2y1 = Cv1v2[ node_l1   + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 || 
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x1   = 0;
          } else {
             var11_x1   = Cv1v1[ node_i_l1 + node_i_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y1   = 0;
          } else {
             var22_y1   = Cv2v2[ node_j_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var12_x1y1 = 0.;  
           } else {
             var12_x1y1 = Cv1v2[ node_i_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var21_x1y1 = 0.;
          }else {
             var21_x1y1 = Cv1v2[ node_i_l2 + node_j_l1 * num_blocks ];
          }
   
          C_vxvx = var11_x2   * x_x1_l1 * x_x1_l2
                 + var11_x12  * x2_x_l1 * x_x1_l2
                 + var11_x21  * x_x1_l1 * x2_x_l2
                 + var11_x1   * x2_x_l1 * x2_x_l2;
          C_vyvy = var22_y2   * y_y1_l1 * y_y1_l2
                 + var22_y12  * y2_y_l1 * y_y1_l2
                 + var22_y21  * y_y1_l1 * y2_y_l2
                 + var22_y1   * y2_y_l1 * y2_y_l2;
          C_vxvy = var12_x2y2 * x_x1_l1 * y_y1_l2
                 + var12_x1y2 * x2_x_l1 * y_y1_l2 
                 + var12_x2y1 * x_x1_l1 * y2_y_l2
                 + var12_x1y1 * x2_x_l1 * y2_y_l2;
          C_vyvx = var21_x2y2 * y_y1_l1 * x_x1_l2
                 + var21_x1y2 * y2_y_l1 * x_x1_l2
                 + var21_x2y1 * y_y1_l1 * x2_x_l2
                 + var21_x1y1 * y2_y_l1 * x2_x_l2;
          Cqq[p1 + p2 * nstrl_ ] =( C_vxvx * costh_l1 * costh_l2
                                  + C_vyvy * sinth_l1 * sinth_l2
                                  + C_vyvx * sinth_l1 * costh_l2
                                  + C_vxvy * costh_l1 * sinth_l2
                                  ) * radius * radius * weight[p1] * weight[p2];
       }
   }

   double Cqq_xi_xi = 0;
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = 0; p1 < nstrl_; ++p1) {
           Cqq_xi_xi += Cqq[p1 + p2 * nstrl_ ];
       }
   }
   /*
   if(unit == FIELD ) {
     double aaa  =     ft_m/day_sec;
     cout << "the variance of Cqq(xi,xi) = " << Cqq_xi_xi/aaa/aaa << ' '<< endl;
   }else {
     cout << "the variance of Cqq(xi,xi) = " << Cqq_xi_xi << ' '<< endl;
   }
   exit(0);
   */
   // At Production Well
   double cell_x_len = prodBlkPtr->getDx();
   double cell_y_len = prodBlkPtr->getDy();
   double cell_z_len = prodBlkPtr->getDz();
   int node, node_left, node_down;
   int iprod = prodBlkPtr->getBlockCtPnt().getI();
   int jprod = prodBlkPtr->getBlockCtPnt().getJ();
   int kprod = prodBlkPtr->getBlockCtPnt().getK();
   double poro = domainPtr->getPoro(iprod, jprod, kprod);
   node      = domainPtr->getGrid()->getIndex(iprod    ,jprod    ,kprod);
   node_left = domainPtr->getGrid()->getIndex(iprod - 1,jprod    ,kprod);
   node_down = domainPtr->getGrid()->getIndex(iprod    ,jprod - 1,kprod);
   
   double Cvx2_vx2 = Cv1v1[ node       + node       * num_blocks ];
   double Cvy2_vy2 = Cv2v2[ node       + node       * num_blocks ];
   double Cvx2_vy2 = Cv1v2[ node       + node       * num_blocks ];
   double Cvx1_vx1, Cvx1_vx2, Cvx2_vx1;
   double Cvy1_vy1, Cvy1_vy2, Cvy2_vy1;
   double Cvx1_vy1, Cvx1_vy2, Cvx2_vy1;
   
   if(iprod == 0 && jprod == 0 ) {
      Cvx1_vy1 = 0.;
   } else {
      Cvx1_vy1 = Cv1v2[ node_left  + node_down  * num_blocks ];
   } 
   if(iprod == 0 ) {
      Cvx1_vx1 = 0.;
      Cvx1_vx2 = 0.;
      Cvx2_vx1 = 0.;
      Cvx1_vy2 = 0.;      
   } else {
      Cvx1_vx1 = Cv1v1[ node_left  + node_left  * num_blocks ];
      Cvx1_vx2 = Cv1v1[ node_left  + node       * num_blocks ];
      Cvx2_vx1 = Cv1v1[ node       + node_left  * num_blocks ]; 
      Cvx1_vy2 = Cv1v2[ node_left  + node       * num_blocks ];      
   }
   if(jprod == 0) {
      Cvy1_vy1 = 0.;
      Cvy1_vy2 = 0.;
      Cvy2_vy1 = 0.;
      Cvx2_vy1 = 0.;
   } else { 
      Cvy1_vy1 = Cv2v2[ node_down  + node_down  * num_blocks ];
      Cvy1_vy2 = Cv2v2[ node_down  + node       * num_blocks ];
      Cvy2_vy1 = Cv2v2[ node       + node_down  * num_blocks ];
      Cvx2_vy1 = Cv1v2[ node       + node_down  * num_blocks ];
   }

   double cqq_total =( ( Cvx2_vx2 + Cvx1_vx1 - Cvx1_vx2 - Cvx2_vx1 ) * cell_y_len * cell_y_len
	              +( Cvy2_vy2 + Cvy1_vy1 - Cvy1_vy2 - Cvy2_vy1 ) * cell_x_len * cell_x_len
		      +( Cvx2_vy2 + Cvx1_vy1 - Cvx1_vy2 - Cvx2_vy1 ) * cell_x_len * cell_y_len * 2.
		   ) * cell_z_len * cell_z_len * poro * poro;
   
   if(unit == FIELD) { 
      double aaa2 = pow(ft_m/day_sec, 2);
      aaa2        = stbPday_m3Psec * stbPday_m3Psec;
      cout <<"variance = " << cqq_total /aaa2 << " by Field Unit" << endl;
      //cout << "node = "      << node      << endl;
      //cout << "node_left = " << node_left << endl;
      //cout << "node_down = " << node_down << endl;
      /*
      cout << "Cvx2_vx2 = " << Cvx2_vx2 /aaa2 << endl; 
      cout << "Cvx1_vx1 = " << Cvx1_vx1 /aaa2 << endl;
      cout << "Cvx2_vx1 = " << Cvx2_vx1 /aaa2 << endl; 
      cout << "Cvx1_vx2 = " << Cvx1_vx2 /aaa2 << endl;
   
      cout << "Cvy2_vy2 = " << Cvy2_vy2 /aaa2 << endl; 
      cout << "Cvy1_vy1 = " << Cvy1_vy1 /aaa2 << endl;
      cout << "Cvy2_vy1 = " << Cvy2_vy1 /aaa2 << endl; 
      cout << "Cvy1_vy2 = " << Cvy1_vy2 /aaa2 << endl;
   
      cout << "Cvx2_vy2 = " << Cvx2_vy2 /aaa2 << endl; 
      cout << "Cvx1_vy1 = " << Cvy1_vy1 /aaa2 << endl;
      cout << "Cvx2_vy1 = " << Cvx2_vy1 /aaa2 << endl; 
      cout << "Cvx1_vy2 = " << Cvx1_vy2 /aaa2 << endl;
      */
      //exit(0);
   } else {
      cout <<"variance = " << cqq_total << endl;
   }
   //exit(0);
   /* 
   // scale back -- not applicable!!!
   double sum = 0;
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = 0; p1 < nstrl_; ++p1) {
	   Cqq[p1 + p2 * nstrl_ ] = Cqq[p1 + p2 * nstrl_ ] * cqq_total/Cqq_xi_xi;
           sum += Cqq[p1 + p2 * nstrl_ ];
       }
   }
   if(unit == FIELD) { 
      double aaa2 = pow(ft_m/day_sec, 2);	   
      cout <<"variance = " << sum / aaa2 << " by Field Unit" << endl; 
   } else {
      cout <<"variance = " << sum << endl;
   }
   */
   //exit(0);

   /*	   
   cout << node << ' ' << node_left << ' ' << node_down << endl; 
   cout << Cvx2_vx2 << ' '
	<< Cvx1_vx1 << ' '
	<< Cvx1_vx2 << ' '
	<< Cvx2_vx1 << endl;
   cout << Cvy2_vy2 << ' '
	<< Cvy1_vy1 << ' '
	<< Cvy1_vy2 << ' '
	<< Cvy2_vy1 << endl;
   cout << Cvx2_vy2 << ' '
	<< Cvx1_vy1 << ' '
	<< Cvx1_vy2 << ' '
	<< Cvx2_vy1 << endl;
	*/
    
   /* //check Symmetry and works!
   cout << setiosflags(ios::fixed | ios::showpoint) << setprecision(6);
   for(int p2 = 0; p2 < nstrl_; ++p2){
       for(int p1 = 0; p1 < nstrl_; ++p1) {
          cout << p1 << ' '
               << p2 << ' ' 
               << Cqq[p1 + p2 * nstrl_ ] << ' '
               << Cqq[p2 + p1 * nstrl_ ] << ' '
               << Cqq[p1 + p2 * nstrl_ ] - Cqq[p2 + p1 * nstrl_ ] << endl;
       }
       cout << endl;
   }
   */
};

void ParticleTrack::calcMoments(Particle *particle1, Particle *particle2,
   double *Vxi_Vxi, double *Vxi_Eta, double *Eta_Eta, 
   double *Vxi_Veta, double *Veta_Veta){

   int len2 = particle2->getLength();
   int len1 = particle1->getLength(); 
   Block *block_l1, *block_l2;
   double C_vxvx, C_vyvy, C_vxvy, C_vyvx;
   int node_l1, node_l2;
   int node_i_l1, node_i_l2;
   int node_j_l1, node_j_l2;
   double var11_x , var11_x1,   var11_x2,   var11_x12,  var11_x21;
   double var22_y , var22_y1,   var22_y2,   var22_y12,  var22_y21;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double var21_xy, var21_x1y1, var21_x1y2, var21_x2y1, var21_x2y2;
   int num_blocks  = domainPtr->getNumBlocks();
   int l1l2;
   double  costh_l2, sinth_l2, costh_l1, sinth_l1;
   double  x_x1_l2, x2_x_l2, y_y1_l2, y2_y_l2;
   double  x_x1_l1, x2_x_l1, y_y1_l1, y2_y_l1;
   for(int l2 = 0; l2 < len2; ++l2) {
        block_l2 = domainPtr->point2block( particle2->getTrajectory(l2) );
         node_l2 = particle2->getNode(l2);
       node_i_l2 = particle2->getNode_i(l2);
       node_j_l2 = particle2->getNode_j(l2);
        costh_l2 = particle2->getCosth(l2);
        sinth_l2 = particle2->getSinth(l2);
         x_x1_l2 = particle2->getX_x1(l2);
         x2_x_l2 = particle2->getX2_x(l2);
         y_y1_l2 = particle2->getY_y1(l2);
         y2_y_l2 = particle2->getY2_y(l2);
       for(int l1 = 0; l1 < len1; ++l1) {
            block_l1 = domainPtr->point2block( particle1->getTrajectory(l1) );
             node_l1 = particle1->getNode(l1);
           node_i_l1 = particle1->getNode_i(l1);
           node_j_l1 = particle1->getNode_j(l1);
            costh_l1 = particle1->getCosth(l1);
            sinth_l1 = particle1->getSinth(l1);
             x_x1_l1 = particle1->getX_x1(l1);
             x2_x_l1 = particle1->getX2_x(l1);
             y_y1_l1 = particle1->getY_y1(l1);
             y2_y_l1 = particle1->getY2_y(l1);
              
          var11_x2   = Cv1v1[ node_l1   + node_l2   * num_blocks ];
          var22_y2   = Cv2v2[ node_l1   + node_l2   * num_blocks ];
          var12_x2y2 = Cv1v2[ node_l1   + node_l2   * num_blocks ];
          var21_x2y2 = Cv1v2[ node_l2   + node_l1   * num_blocks ];
          if(block_l1->getBlockLbPnt().getI() == 0 ) {
             var11_x12  = var12_x1y2 = 0;
          } else { 
             var11_x12  = Cv1v1[ node_i_l1 + node_l2  * num_blocks ];
             var12_x1y2 = Cv1v2[ node_i_l1 + node_l2  * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x21  = var21_x2y1 = 0;
          } else {
             var11_x21  = Cv1v1[ node_l1   + node_i_l2 * num_blocks ];
             var21_x2y1 = Cv1v2[ node_i_l2 + node_l1   * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ) {
             var22_y12  = var21_x1y2 = 0;          
          } else {
             var22_y12  = Cv2v2[ node_j_l1 + node_l2   * num_blocks ];
             var21_x1y2 = Cv1v2[ node_l2   + node_j_l1 * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y21  = var12_x2y1 = 0;          
          } else {
             var22_y21  = Cv2v2[ node_l1   + node_j_l2 * num_blocks ];
             var12_x2y1 = Cv1v2[ node_l1   + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 || 
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x1   = 0;
          } else {
             var11_x1   = Cv1v1[ node_i_l1 + node_i_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y1   = 0;
          } else {
             var22_y1   = Cv2v2[ node_j_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var12_x1y1 = 0.;  
           } else {
             var12_x1y1 = Cv1v2[ node_i_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var21_x1y1 = 0.;
          }else {
             var21_x1y1 = Cv1v2[ node_i_l2 + node_j_l1 * num_blocks ];
          }
   
          C_vxvx = var11_x2   * x_x1_l1 * x_x1_l2
                 + var11_x12  * x2_x_l1 * x_x1_l2
                 + var11_x21  * x_x1_l1 * x2_x_l2
                 + var11_x1   * x2_x_l1 * x2_x_l2;
          C_vyvy = var22_y2   * y_y1_l1 * y_y1_l2
                 + var22_y12  * y2_y_l1 * y_y1_l2
                 + var22_y21  * y_y1_l1 * y2_y_l2
                 + var22_y1   * y2_y_l1 * y2_y_l2;
          C_vxvy = var12_x2y2 * x_x1_l1 * y_y1_l2
                 + var12_x1y2 * x2_x_l1 * y_y1_l2 
                 + var12_x2y1 * x_x1_l1 * y2_y_l2
                 + var12_x1y1 * x2_x_l1 * y2_y_l2;
          C_vyvx = var21_x2y2 * y_y1_l1 * x_x1_l2
                 + var21_x1y2 * y2_y_l1 * x_x1_l2
                 + var21_x2y1 * y_y1_l1 * x2_x_l2
                 + var21_x1y1 * y2_y_l1 * x2_x_l2;
          l1l2 = l1 + l2 * len1;
          Vxi_Vxi[ l1l2 ]   = C_vxvx * costh_l1 * costh_l2
                            + C_vyvx * sinth_l1 * costh_l2
                            + C_vxvy * costh_l1 * sinth_l2
                            + C_vyvy * sinth_l1 * sinth_l2;
          Veta_Veta[ l1l2 ] = C_vxvx * sinth_l1 * sinth_l2
                            - C_vyvx * costh_l1 * sinth_l2
                            - C_vxvy * sinth_l1 * costh_l2
                            + C_vyvy * costh_l1 * costh_l2;
          Vxi_Veta[ l1l2 ] =- C_vxvx * costh_l1 * sinth_l2
                            - C_vyvx * sinth_l1 * sinth_l2
                            + C_vxvy * costh_l1 * costh_l2
                            + C_vyvy * sinth_l1 * costh_l2;
       }
   }

   //(2) Vxi_Eta
   int l1l2m1;
   for(int l1 = 0; l1 < len1; ++l1) {
       l1l2 = l1;
       Vxi_Eta[l1l2] = 0.;
   }
   
   for(int l2 = 1; l2 < len2; ++l2) {
       for(int l1 = 0; l1 < len1; ++l1) {
          l1l2   = l1 +  l2      * len1;
          l1l2m1 = l1 + (l2 - 1) * len1;
          Vxi_Eta[l1l2] = Vxi_Eta[l1l2m1] 
                       +( Vxi_Veta[l1l2]  /particle2->getVAvg(l2) 
                       +  Vxi_Veta[l1l2m1]/particle2->getVAvg(l2-1)
                        ) /2. * particle2->getDs(l2);
       }
   }

   //(3) eta_eta
   int l1m1l2, l1m1l2m1;
   for(int l2 = 0; l2 < len2; ++l2) {
       l1l2 = l2 * len1;
       Eta_Eta[l1l2] = 0.;
   }
   for(int l1 = 0; l1 < len1; ++l1) {
       l1l2     = l1;
       Eta_Eta[l1l2] = 0.;
   }
   double tmp;
   for(int l2 = 1; l2 < len2; ++l2) {
       for(int l1 = 1; l1 < len1; ++l1) {
         l1l2     = l1     + l2 * len1;   
         l1m1l2   = l1 - 1 + l2 * len1;
         l1l2m1   = l1     + (l2 - 1) * len1;   
         l1m1l2m1 = l1 - 1 + (l2 - 1) * len1;
         tmp = ( Veta_Veta[l1l2]   /particle1->getVAvg(l1)
                                   /particle2->getVAvg(l2)   
              + Veta_Veta[l1m1l2]  /particle1->getVAvg(l1-1)
                                   /particle2->getVAvg(l2)
              + Veta_Veta[l1l2m1]  /particle1->getVAvg(l1)
                                   /particle2->getVAvg(l2-1)
              + Veta_Veta[l1m1l2m1]/particle1->getVAvg(l1-1)
                                   /particle2->getVAvg(l2-1)
              )/4.*particle1->getDs(l1)*particle2->getDs(l2);
          Eta_Eta[l1l2] = Eta_Eta[l1m1l2] + Eta_Eta[l1l2m1] 
                        + tmp -  Eta_Eta[l1m1l2m1];
       }
   }       
}

void ParticleTrack::calcMoments(Particle *particle1, Particle *particle2,
   double *Vxi_Veta, double *Vxi_Eta){ 

   int len2 = particle2->getLength();
   int len1 = particle1->getLength(); 
   Block *block_l1, *block_l2;
   double C_vxvx, C_vyvy, C_vxvy, C_vyvx;
   int node_l1, node_l2;
   int node_i_l1, node_i_l2;
   int node_j_l1, node_j_l2;
   double var11_x , var11_x1,   var11_x2,   var11_x12,  var11_x21;
   double var22_y , var22_y1,   var22_y2,   var22_y12,  var22_y21;
   double var12_xy, var12_x1y1, var12_x1y2, var12_x2y1, var12_x2y2;
   double var21_xy, var21_x1y1, var21_x1y2, var21_x2y1, var21_x2y2;
   int num_blocks  = domainPtr->getNumBlocks();
   int l1l2;
   double  costh_l2, sinth_l2, costh_l1, sinth_l1;
   double  x_x1_l2, x2_x_l2, y_y1_l2, y2_y_l2;
   double  x_x1_l1, x2_x_l1, y_y1_l1, y2_y_l1;
   for(int l2 = 0; l2 < len2; ++l2) {
        block_l2 = domainPtr->point2block( particle2->getTrajectory(l2) );
         node_l2 = particle2->getNode(l2);
       node_i_l2 = particle2->getNode_i(l2);
       node_j_l2 = particle2->getNode_j(l2);
        costh_l2 = particle2->getCosth(l2);
        sinth_l2 = particle2->getSinth(l2);
         x_x1_l2 = particle2->getX_x1(l2);
         x2_x_l2 = particle2->getX2_x(l2);
         y_y1_l2 = particle2->getY_y1(l2);
         y2_y_l2 = particle2->getY2_y(l2);
       for(int l1 = 0; l1 < len1; ++l1) {
            block_l1 = domainPtr->point2block( particle1->getTrajectory(l1) );
             node_l1 = particle1->getNode(l1);
           node_i_l1 = particle1->getNode_i(l1);
           node_j_l1 = particle1->getNode_j(l1);
            costh_l1 = particle1->getCosth(l1);
            sinth_l1 = particle1->getSinth(l1);
             x_x1_l1 = particle1->getX_x1(l1);
             x2_x_l1 = particle1->getX2_x(l1);
             y_y1_l1 = particle1->getY_y1(l1);
             y2_y_l1 = particle1->getY2_y(l1);
              
          var11_x2   = Cv1v1[ node_l1   + node_l2   * num_blocks ];
          var22_y2   = Cv2v2[ node_l1   + node_l2   * num_blocks ];
          var12_x2y2 = Cv1v2[ node_l1   + node_l2   * num_blocks ];
          var21_x2y2 = Cv1v2[ node_l2   + node_l1   * num_blocks ];
          if(block_l1->getBlockLbPnt().getI() == 0 ) {
             var11_x12  = var12_x1y2 = 0;
          } else { 
             var11_x12  = Cv1v1[ node_i_l1 + node_l2  * num_blocks ];
             var12_x1y2 = Cv1v2[ node_i_l1 + node_l2  * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x21  = var21_x2y1 = 0;
          } else {
             var11_x21  = Cv1v1[ node_l1   + node_i_l2 * num_blocks ];
             var21_x2y1 = Cv1v2[ node_i_l2 + node_l1   * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ) {
             var22_y12  = var21_x1y2 = 0;          
          } else {
             var22_y12  = Cv2v2[ node_j_l1 + node_l2   * num_blocks ];
             var21_x1y2 = Cv1v2[ node_l2   + node_j_l1 * num_blocks ];
          }
          if(block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y21  = var12_x2y1 = 0;          
          } else {
             var22_y21  = Cv2v2[ node_l1   + node_j_l2 * num_blocks ];
             var12_x2y1 = Cv1v2[ node_l1   + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 || 
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var11_x1   = 0;
          } else {
             var11_x1   = Cv1v1[ node_i_l1 + node_i_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var22_y1   = 0;
          } else {
             var22_y1   = Cv2v2[ node_j_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getI() == 0 ||
             block_l2->getBlockLbPnt().getJ() == 0 ) {
             var12_x1y1 = 0.;  
           } else {
             var12_x1y1 = Cv1v2[ node_i_l1 + node_j_l2 * num_blocks ];
          }
          if(block_l1->getBlockLbPnt().getJ() == 0 ||
             block_l2->getBlockLbPnt().getI() == 0 ) {
             var21_x1y1 = 0.;
          }else {
             var21_x1y1 = Cv1v2[ node_i_l2 + node_j_l1 * num_blocks ];
          }
   
          C_vxvx = var11_x2   * x_x1_l1 * x_x1_l2
                 + var11_x12  * x2_x_l1 * x_x1_l2
                 + var11_x21  * x_x1_l1 * x2_x_l2
                 + var11_x1   * x2_x_l1 * x2_x_l2;
          C_vyvy = var22_y2   * y_y1_l1 * y_y1_l2
                 + var22_y12  * y2_y_l1 * y_y1_l2
                 + var22_y21  * y_y1_l1 * y2_y_l2
                 + var22_y1   * y2_y_l1 * y2_y_l2;
          C_vxvy = var12_x2y2 * x_x1_l1 * y_y1_l2
                 + var12_x1y2 * x2_x_l1 * y_y1_l2 
                 + var12_x2y1 * x_x1_l1 * y2_y_l2
                 + var12_x1y1 * x2_x_l1 * y2_y_l2;
          C_vyvx = var21_x2y2 * y_y1_l1 * x_x1_l2
                 + var21_x1y2 * y2_y_l1 * x_x1_l2
                 + var21_x2y1 * y_y1_l1 * x2_x_l2
                 + var21_x1y1 * y2_y_l1 * x2_x_l2;
          l1l2 = l1 + l2 * len1;
          Vxi_Veta[ l1l2 ] =- C_vxvx * costh_l1 * sinth_l2
                            - C_vyvx * sinth_l1 * sinth_l2
                            + C_vxvy * costh_l1 * costh_l2
                            + C_vyvy * sinth_l1 * costh_l2;
       }
   }

   //(2) Vxi_Eta
   int l1l2m1;
   for(int l1 = 0; l1 < len1; ++l1) {
       l1l2 = l1;
       Vxi_Eta[l1l2] = 0.;
   }
   
   for(int l2 = 1; l2 < len2; ++l2) {
       for(int l1 = 0; l1 < len1; ++l1) {
          l1l2   = l1 +  l2      * len1;
          l1l2m1 = l1 + (l2 - 1) * len1;
          Vxi_Eta[l1l2] = Vxi_Eta[l1l2m1] 
                       +( Vxi_Veta[l1l2]  /particle2->getVAvg(l2) 
                       +  Vxi_Veta[l1l2m1]/particle2->getVAvg(l2-1)
                        ) /2. * particle2->getDs(l2);
       }
   }
}

void ParticleTrack::printTravelTime() {
   ofstream fileOu("TravelTime.out", ios::out);
   if (!fileOu){
       cerr << "\n File does not exist !\n";
       exit(EXIT_FAILURE);
   }
   for(int i = 0; i < nstrl_; ++i){
       fileOu<<particles[i]->getLastTravelTime() << ' '; 
       fileOu<<particles[i]->getLastTravelTVar() << endl;
   }
   fileOu.close();
}

void ParticleTrack::printLogTravelTime() {
   ofstream fileOu("LogTravelTime.out", ios::out);
   if (!fileOu){
       cerr << "\n File does not exist !\n";
       exit(EXIT_FAILURE);
   }
   for(int i = 0; i < nstrl_; ++i){
       fileOu << trvl_avg_log[i]<<' ' 
              << trvl_var_log[i]<<' '
	      << sqrt( trvl_var_log[i]) << endl;
   }

   fileOu.close();
}

