/*
 * File: PVar_Sensitivity.cc
 * ----------------------------------
 * Implementation for PVar_Sensitivity class
 */

#include "PVar_Sensitivity.h"

//===== PVar_Sensitivity( const P0Eqn &P0ee ) =====
PVar_Sensitivity::PVar_Sensitivity(const P0Eqn &P0ee, CYPEqn *CYPee, CPPEqn *CPPee,
                  MasterPoint *mpts, PressMeasure *prsm, KrigWeigh *krig)
: Eqn( P0ee ), P0e( &P0ee ), CYPe(CYPee), CPPe(CPPee),
  mptsPtr(mpts), prsmPtr(prsm), krigPtr(krig) {
 
  debug = false; 
  if(debug) 
     cout << "PVar_Sensitivity::PVar_Sensitivity( )"<<endl;
 
  bigNumber = 1.0e14;  
  
  num_meas = prsmPtr->getLength();

  p0  = P0e->getP0();  
  
  pVar_Sensitivity = new double[ num_meas * mptsPtr->getLength() ];
  for(int i = 0; i < num_meas * mptsPtr->getLength(); ++i) {
      pVar_Sensitivity[i] = 0.0;
  }
  
  rho_YP      = new double[nNode];
  std_PP      = new double[nNode]; 
  std_YY      = new double[nNodePerm];
  std_PP_UnCd = new double[nNode]; 
  std_YY_UnCd = new double[nNodePerm];
  setStdP();
}

//===== ~PVar_Sensitivity() =====
PVar_Sensitivity::~PVar_Sensitivity() {
  if(debug) 
     cout << "PVar_Sensitivity::~PVar_Sensitivity( )"<<endl;
  delete[] pVar_Sensitivity;
  delete[] rho_YP;
  delete[] std_PP;
  delete[] std_YY;
  delete[] std_PP_UnCd;
  delete[] std_YY_UnCd;  
}

void PVar_Sensitivity::setStdP() {
  for(int i = 0; i < nNode; ++i) {
      std_PP[i] = CPPe->getUnCdPStd(i);
      std_PP_UnCd[i] = std_PP[i];
  }
  double sum = 0;
  for(int kc = 0; kc < nzPerm; kc++) {
      for(int jc = 0; jc < nyPerm; jc++) {
          for(int ic = 0; ic < nxPerm; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
              std_YY[ijkc] = sqrt( permPtr->getYVar(ic,jc,kc) );
              std_YY_UnCd[ijkc] = std_YY[ijkc];
              sum += std_YY[ijkc];
          }
      }
  }
  sum /= (nzPerm * nyPerm * nxPerm);
  if(debug)
     cout << "the global std of unconditonal Y = " 
          << sum << endl;
        
  // IMPORTANT NOTE, rho_YP is only valid at ijkc = 0 here !!!
  rho_k   = 0;
  rho_j   = 0;
  rho_i   = 0;
  rho_ijk = g.getPermIndex(rho_i, rho_j, rho_k);
  for(int i = 0; i < nNode; ++i) {
      if(fabs(std_PP[i]) > 0.00001) {
         rho_YP[i] = CYPe->getCYP(i,rho_ijk)
                     / std_PP_UnCd[i]/ std_YY_UnCd[rho_ijk];
      }
      else {
         rho_YP[i] = 0.0;
      }
  } 
}

void PVar_Sensitivity::calcSensitity(int nWells, const Well *w, 
                       const Control &c, Solver &s) {
   	
   for(int mst = 0; mst < mptsPtr->getLength(); mst++) {
       cout << mptsPtr->getIm(mst) << ' '
	    << mptsPtr->getJm(mst) << ' '
	    << mptsPtr->getKm(mst) << ' '
	    << mptsPtr->getIJK(mst) << ' '   
	    << endl;
       ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
       CYPe->solve(nWells, w, c, s, 0);
       CPPe->solve(nWells, w, c, s, 0);
   }
   exit(0);

   
   for(int i = 0; i < nNode; ++i) {
      std_PP[i] = CPPe->getUnCdPStd(i);
   }
}

//=====  calcSensitity() =====
void PVar_Sensitivity::calcSensitity(int nWells, const Well *w, 
                       const Control &c, Solver &s, int iter ) {
                             
  // --- point to P0's derivatives for current timestep -------
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, 0) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, 0) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, 0) ;

  for(int kc = 0; kc < nzPerm; kc++) {
      for(int jc = 0; jc < nyPerm; jc++) {
          for(int ic = 0; ic < nxPerm; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
              std_YY[ijkc] = sqrt( permPtr->getYVar(ic,jc,kc) );
          }
      }
  }

  //ofstream os("PVar_CPP.out", ios::out);
  int calcType = 2; 
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      int ic = mptsPtr->getIm( m );
      int jc = mptsPtr->getJm( m );
      int kc = mptsPtr->getKm( m );
      int ijkc = g.getPermIndex(ic, jc, kc);
      if(debug) 
         cout << ic <<' '<< jc <<' '<< kc <<' '<<m <<endl;
      for(int i = 0; i < nNode; ++i) {
          if(fabs(std_PP_UnCd[i]) > 0.00001) {
             rho_YP[i] = CYPe->getCYP(i,ijkc)
                         / std_PP_UnCd[i]/ std_YY_UnCd[ijkc];
          }
          else {
             rho_YP[i] = 0.0;
          }
      } 
      calcARHO(c.getBType(), c.getBCond(), nWells, w, calcType);
      /*
      for(int i = 0; i < nNode; i++) os << diagA[i] << endl;
      os << endl;
      for(int i = 0; i < nNode - 1; i++) 
          os<< offDAXp[i]<< ' '<< offDAXn[i]<<endl;
      os << endl;
      for(int i = 0; i < nNode - nx; i++) 
          os<< offDAYp[i]<< ' '<< offDAYn[i]<<endl ;
      exit(0);
      */
      calcSensiRHS(m, ic, jc, kc, c.getBType(), c.getBCond());
      //calcRHS(ic, jc, kc, c.getBType(), c.getBCond());
      /*
      for(int i = 0; i < nNode; ++i) {
          os << RHS[i] <<endl;
          if( ((i+1)%nx) ==0 ) os << endl;
      }
      exit(0);
      */
      s.solve( true, *this, 1, RHS );
      /*
      for(int i = 0; i < nNode; ++i) {
          os << RHS[i] << ' ' << std_PP[i]<<endl;
          if( ((i+1)%nx) ==0 ) os << endl;
      }
      */
      for(int i = 0; i < num_meas; ++i) {
          int ijk = g.getIndex( prsmPtr->getIpmeas(i),
                                prsmPtr->getJpmeas(i),
                                prsmPtr->getKpmeas(i)
                              );
          pVar_Sensitivity[ i + m * num_meas] = RHS[ ijk ];
      }
      //exit(0);
  } 

  if(debug) { output();
     cout << "Done at PVar_Sensitivity::solve () " <<endl;
  }
}

void PVar_Sensitivity::calcCondPStdByFullEqn(int nWells, const Well *w, 
                const Control &c, Solver &s, int iter) {
   ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
   CYPe->solve(nWells, w, c, s, 0);
   CPPe->solve(nWells, w, c, s, 0);
   for(int i = 0; i < nNode; ++i) {
      std_PP[i] = CPPe->getUnCdPStd(i);
   }
}

void PVar_Sensitivity::calcSensitivityByFullEqn(int nWells, const Well *w, 
   const Control &c, Solver &s, int iter) {

   int i_cond = 2;
   int j_cond = 2;
   int k_cond = 0;
   int i2_cond = g.toIPerm( i_cond );
   int j2_cond = g.toJPerm( j_cond );
   int k2_cond = g.toKPerm( k_cond );
   int ijk2_left  = g.getPermIndex(i2_cond - 1, j2_cond, k2_cond);
   int ijk2_right = g.getPermIndex(i2_cond + 1, j2_cond, k2_cond);
   int ijk2_down  = g.getPermIndex(i2_cond    , j2_cond - 1, k2_cond);
   int ijk2_upper = g.getPermIndex(i2_cond    , j2_cond + 1, k2_cond);

   int ijk_ltu = g.getIndex(i_cond - 1, j_cond + 1, k_cond);   
   int ijk_u   = g.getIndex(i_cond    , j_cond + 1, k_cond);
   int ijk_rtu = g.getIndex(i_cond + 1, j_cond + 1, k_cond);   
   int ijk_lt  = g.getIndex(i_cond - 1, j_cond    , k_cond);
   int ijk     = g.getIndex(i_cond    , j_cond    , k_cond);
   int ijk_rt  = g.getIndex(i_cond + 1, j_cond    , k_cond);
   int ijk_ltd = g.getIndex(i_cond - 1, j_cond - 1, k_cond);
   int ijk_d   = g.getIndex(i_cond    , j_cond - 1, k_cond);
   int ijk_rtd = g.getIndex(i_cond + 1, j_cond - 1, k_cond);
   
   permPtr->setKAvg(ijk2_left );
   permPtr->setKAvg(ijk2_right);
   permPtr->setKAvg(ijk2_down );
   permPtr->setKAvg(ijk2_upper);
   ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
   CYPe->solve(nWells, w, c, s, 0);
   CPPe->solve(nWells, w, c, s, 0);
   
   double pvar_old = CPPe->getCPP(ijk, ijk);
   double pvar_new = CPPe->getCPP(ijk, ijk);

   double dpvar_x  = CPPe->getCPP(ijk_lt , ijk_lt ) - pvar_old;
   double dpvar_y  = CPPe->getCPP(ijk_d  , ijk_d  ) - pvar_old;
   double dpvar_xy = CPPe->getCPP(ijk_ltd, ijk_ltd) - pvar_old;
	          
   int    i_chosen = 0, i_chosen_2nd = 0;
   bool   found, found_2nd;
   double dperm      = 0.1, dpres;
   double dperm_2nd  = 0.1, dpres_2nd;
   double perm_value = 1.0;
   double diff;
   cout<< " Unconditional : " << endl;
   cout<< sqrt( CPPe->getCPP(ijk_ltd, ijk_ltd) ) << ' '
       << sqrt( CPPe->getCPP(ijk    , ijk    ) ) << ' '
       << sqrt( CPPe->getCPP(ijk_rtu, ijk_rtu) ) << ' '
       << endl;
   
   for(int iter = 0; iter < 50000; ++iter) {
       found     = false;
       found_2nd = false;
       diff = 0;
       for(int i = 0; i < nNodePerm; ++i) {
           double perm_value_old = permPtr->getYVar(i);
           double perm_value_new = perm_value_old  - dperm;
           if( perm_value_new > 1.) perm_value_new = 1.;
	   if( perm_value_new < 0.) perm_value_new = 0.;

	   permPtr->setYVar(i, perm_value_new);
	   if( i == ijk2_left || i == ijk2_right ||
	       i == ijk2_down || i == ijk2_upper ) {
	       permPtr->setKAvg(ijk2_left );
               permPtr->setKAvg(ijk2_right);
               permPtr->setKAvg(ijk2_down );
               permPtr->setKAvg(ijk2_upper);
	       ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
	   }

           CYPe->solve(nWells, w, c, s, 0);
           CPPe->solve(nWells, w, c, s, 0);
	   
	   if( (CPPe->getCPP(ijk_ltd, ijk_ltd) - CPPe->getCPP(ijk, ijk)) > dpvar_xy ) {
              i_chosen   = i;
	      found      = true;
	      pvar_new   = CPPe->getCPP(ijk, ijk);
	      dpvar_xy   = CPPe->getCPP(ijk_ltd, ijk_ltd) - pvar_new;
	      dpres      = pvar_old - pvar_new;
              dperm      = perm_value_old - perm_value_new;
           } 
	   /*
	   else if
             (CPPe->getCPP(ijk_lt, ijk_lt) - CPPe->getCPP(ijk, ijk) > dpvar_x &&
	      CPPe->getCPP(ijk_d , ijk_d ) - CPPe->getCPP(ijk, ijk) > dpvar_y ) {
	      i_chosen   = i;
	      found      = true;
	      pvar_new   = CPPe->getCPP(ijk, ijk);
	      dpvar_x    = CPPe->getCPP(ijk_lt, ijk_lt) - pvar_new;
	      dpvar_y    = CPPe->getCPP(ijk_d , ijk_d ) - pvar_new;
	      dpres      = pvar_old - pvar_new;
              dperm      = perm_value_old - perm_value_new;
	   } */
	   
	   if(pvar_old - CPPe->getCPP(ijk, ijk) > diff ) {
	      i_chosen_2nd = i;
	      found_2nd    = true;
	      diff         = pvar_old - CPPe->getCPP(ijk, ijk);
	      pvar_new     = CPPe->getCPP(ijk, ijk);
	      dpres_2nd    = pvar_old - pvar_new;
              dperm_2nd    = perm_value_old - perm_value_new;
	   }
	   
	   permPtr->setYVar(i, perm_value_old);
       }
       if( found ) {
          perm_value = permPtr->getYVar(i_chosen) - dperm/dpres * pvar_old;
          if(perm_value < 0) perm_value = 0.;
	  if(perm_value > 1) perm_value = 1.;
	  cout << "Iter = "  << iter <<' ' 
               << i_chosen   << ' ' 
	       << perm_value << endl;
          permPtr->setYVar(i_chosen, perm_value);
	  if( i_chosen == ijk2_left || i_chosen == ijk2_right ||
	      i_chosen == ijk2_down || i_chosen == ijk2_upper 
	    ) {
	      permPtr->setKAvg(ijk2_left );
              permPtr->setKAvg(ijk2_right);
              permPtr->setKAvg(ijk2_down );
              permPtr->setKAvg(ijk2_upper);
	      ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
	  }
	  
          CYPe->solve(nWells, w, c, s, 0);
          CPPe->solve(nWells, w, c, s, 0);

          cout << sqrt( CPPe->getCPP(20, 20) ) << ' '
	       << sqrt( CPPe->getCPP(21, 21) ) << ' '
               << sqrt( CPPe->getCPP(22, 22) ) << ' '
               << sqrt( CPPe->getCPP(23, 23) ) << ' '
	       << sqrt( CPPe->getCPP(24, 24) ) <<endl;
          cout << sqrt( CPPe->getCPP(15, 15) ) << ' '
               << sqrt( CPPe->getCPP(16, 16) ) << ' '
	       << sqrt( CPPe->getCPP(17, 17) ) << ' '
	       << sqrt( CPPe->getCPP(18, 18) ) << ' '
	       << sqrt( CPPe->getCPP(19, 19) ) <<endl;
          cout << sqrt( CPPe->getCPP(10, 10) ) << ' '
	       << sqrt( CPPe->getCPP(11, 11) ) << ' '
               << sqrt( CPPe->getCPP(12, 12) ) << ' '
               << sqrt( CPPe->getCPP(13, 13) ) << ' '
	       << sqrt( CPPe->getCPP(14, 14) ) <<endl;
          cout << sqrt( CPPe->getCPP( 5,  5) ) << ' '
               << sqrt( CPPe->getCPP( 6,  6) ) << ' '
               << sqrt( CPPe->getCPP( 7,  7) ) << ' '
               << sqrt( CPPe->getCPP( 8,  8) ) << ' '
	       << sqrt( CPPe->getCPP( 9,  9) ) << endl;
          cout << sqrt( CPPe->getCPP( 0,  0) ) << ' '
               << sqrt( CPPe->getCPP( 1,  1) ) << ' '
               << sqrt( CPPe->getCPP( 2,  2) ) << ' '
               << sqrt( CPPe->getCPP( 3,  3) ) << ' '
               << sqrt( CPPe->getCPP( 4,  4) ) << endl;
	  
       }  
       else if( found_2nd) {
          perm_value = permPtr->getYVar(i_chosen_2nd) - dperm_2nd/dpres_2nd * pvar_old;
          if(perm_value < 0) perm_value = 0.;
	  if(perm_value > 1) perm_value = 1.;
	  cout << "Iter_2nd = "  << iter <<' ' 
               << i_chosen_2nd        << ' ' 
	       << dperm_2nd/dpres_2nd << ' '
	       << perm_value   << endl;
          permPtr->setYVar(i_chosen_2nd, perm_value);
	  if( i_chosen_2nd == ijk2_left || i_chosen_2nd == ijk2_right ||
	      i_chosen_2nd == ijk2_down || i_chosen_2nd == ijk2_upper 
	    ) {
	      permPtr->setKAvg(ijk2_left );
              permPtr->setKAvg(ijk2_right);
              permPtr->setKAvg(ijk2_down );
              permPtr->setKAvg(ijk2_upper);
	      ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
	  }
          CYPe->solve(nWells, w, c, s, 0);
          CPPe->solve(nWells, w, c, s, 0);
          cout << sqrt( CPPe->getCPP(16, 16) )<< ' '
	       << sqrt( CPPe->getCPP(17, 17) )<< ' '
	       << sqrt( CPPe->getCPP(18, 18) )<< endl;
          cout << sqrt( CPPe->getCPP(11, 11) )<< ' '
               << sqrt( CPPe->getCPP(12, 12) )<< ' '
               << sqrt( CPPe->getCPP(13, 13) )<< endl;
          cout << sqrt( CPPe->getCPP( 6,  6) )<< ' '
               << sqrt( CPPe->getCPP( 7,  7) )<< ' '
               << sqrt( CPPe->getCPP( 8,  8) )<< endl;
       } 
       else {
          exit(0);
       }   
   }
}

void PVar_Sensitivity::calcCondPStd(int nWells, const Well *w, 
                       const Control &c, Solver &s, double* solu) {
  
  //( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
      
  // --- point to P0's derivatives for current timestep -----
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, 0) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, 0) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, 0) ;
  
  // failed try
  // calcSensitivityByFullEqn(nWells, w, c, s, 0);
  // exit(0);
  
  int kc = 0;
  int jc = 0;
  int ic = 0;
  int i_rho = 1;
  if(i_rho == 1) {
     for(int i = 0; i < nNode; ++i) 
         rho_YP[i] = CPPe->getUnCdRho(i);
         //rho_YP[12] = 1.;
  } else if(i_rho == 2)   {
     for(int i = 0; i < nNode; ++i) 
         rho_YP[i] = CPPe->getCondRho(i);
     rho_YP[12] = 1.; 
     //diagA[12] = 1.0e14;
  } else {
     cout << "no implementation!" << endl;
     exit(0);
  }

  calcARHO(c.getBType(), c.getBCond(), nWells, w, 2);
  int choice = 3;
  if( choice == 1) {
    calcRHS( ic, jc, kc, c.getBType(), c.getBCond() );
  } else if(choice == 2) {
    calcConstrain(ic, jc, kc, solu);
    calcRHS( ic, jc, kc, c.getBType(), c.getBCond() );
  } else if(choice == 3) {
    double *ainv = new double [nNode * nNode];
    double *am1  = new double [nNode];
    double *am2  = new double [nNode];
    double *am3  = new double [nNode];
    double *am4  = new double [nNode];
    double *amR  = new double [nNodePerm];
    double *amO  = new double [nNodePerm];
    double *xSolu= new double [nNodePerm];

    for(int i = 0; i < nNodePerm; ++i) {
	amR[i] = amO[i] = 0.;
    }
    
    int node_cond = 12;
    calcAinv(s, ainv);
    calcAMs(ic, jc, kc,  am1, am2, am3, am4);
    
    //Note that j is the row number !!
    /*
    for(int irow = 12; irow < 13; ++irow) {
        for(int jcol = 0; jcol < nNode; ++jcol) {
           am1[jcol] *= ainv[irow + jcol * nNode];  
           am2[jcol] *= ainv[irow + jcol * nNode];
           am3[jcol] *= ainv[irow + jcol * nNode]; 
           am4[jcol] *= ainv[irow + jcol * nNode];
	}
    }*/
    for(int irow = 1; irow < nNode; ++irow) {
    for(int k = 0; k < nz; k++ ) {
        int k2 = g.toKPerm( k );
        for(int j = 0; j < ny; j++ ) {
            int j2 = g.toJPerm( j );
            int j2_p1 = j2 + 1;
            int j2_n1 = j2 - 1;
            for(int i = 0; i < nx; i++ ) {
                int i2     = g.toIPerm( i );
                int i2_p1  = i2 + 1;
                int i2_n1  = i2 - 1;
                int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
                int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
                int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
                int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
		int ijk      = g.getIndex( i  , j  , k );
		if(irow == node_cond) {
		   amR[ijk2_ip1] += am2[ijk] * ainv[irow + ijk * nNode];
		   amR[ijk2_in1] += am1[ijk] * ainv[irow + ijk * nNode];
		   amR[ijk2_jp1] += am4[ijk] * ainv[irow + ijk * nNode];
		   amR[ijk2_jn1] += am3[ijk] * ainv[irow + ijk * nNode];
		}  else {
		   amO[ijk2_ip1] += am2[ijk] * ainv[irow + ijk * nNode];
		   amO[ijk2_in1] += am1[ijk] * ainv[irow + ijk * nNode];
		   amO[ijk2_jp1] += am4[ijk] * ainv[irow + ijk * nNode];
		   amO[ijk2_jn1] += am3[ijk] * ainv[irow + ijk * nNode];
		}
            }
        }
    }
    }
    
    for(int i = 0; i < nNodePerm; ++i) {
        amO[i] = 1.0;	 
    }

    calcPermStdBySimplx(amR, amO, xSolu);
    for(int i = 0; i < nNodePerm; ++i) {
        permPtr->setYVar(i, xSolu[i]); 
    }
    calcRHS( ic, jc, kc, c.getBType(), c.getBCond() );    
    /* 
    double sum = 0.;
    for(int i = 0; i < nNode; ++i) {
        sum +=  am1[i] + am2[i] + am3[i] + am4[i] ;
    }
    cout << "sum = " << sum << endl;
    sum = 0.;
    double neg = 0., pos = 0;
    for(int i = 0; i < nNodePerm; ++i) {
	if(amR[i] < 0) neg += amR[i];
        if(amR[i] > 0) pos += amR[i];	
	sum += amR[i] ;
    }
    cout << "sum = " << sum << endl;
    cout << "neg = " << neg << endl; 
    cout << "pos = " << pos << endl;
    */
    //exit(0);
    
    delete [] am1;
    delete [] am2;
    delete [] am3;
    delete [] am4;
    delete [] amR;
    delete [] amO;
    delete [] xSolu;
    delete [] ainv;
  }
  /*
  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
  for(int n =  0; n < nNode; ++n) 
      cout << "node = " << n <<' '<< RHS[n] << endl;
  */
  
  // 1 fast solution 2, full solution
  int solu_choice = 2; 
  if(solu_choice == 1) {
     s.solve( true, *this, 1, RHS );
     for(int i = 0; i < nNode; ++i) {
         if(fabs(RHS[i]) > 0.0000001)
            std_PP[i] = RHS[i];
         else 
            std_PP[i] = 0.0; 
     }
  } else if( solu_choice == 2){
     calcCondPStdByFullEqn(nWells, w, c, s, 0);
     for(int i = 0; i < nNode; ++i) {
        if( CPPe->getCPP(i, i) > 0) {
	   std_PP[i] = sqrt(CPPe->getCPP(i, i));
	} else {
	   std_PP[i] = 0;
	}
     }
  }

  cout << "std_PP is updated in Pvar" << endl;

  //simplx_driver1();
  //simplx_driver2();
}

void PVar_Sensitivity::calcAinv(Solver &s, double * ainv) {
    for(int i = 0; i < nNode; ++i) {
	for(int j = 0; j < nNode; ++j) RHS[j] = 0.;
       	RHS[i] = 1.0;
        s.solve( true, *this, 1, RHS );
	for(int j = 0; j < nNode; ++j) {
	    ainv[j + i * nNode] = RHS[j];
	}
    }
}

void PVar_Sensitivity::calcAMs(int i1_perm, int j1_perm, 
  int k1_perm, double *am1, double *am2, double *am3, double *am4) {

  for(int i = 0; i < nNode; ++i) am1[i] = am2[i] = am3[i] = am4[i] = 0.;
  double permx_interface[3];
  double permy_interface[3];
  for(int k = 0; k < nz; k++ ) {
      int k2 = g.toKPerm( k );
      for(int j = 0; j < ny; j++ ) {
          int j2 = g.toJPerm( j );
          int j2_p1 = j2 + 1;
          int j2_n1 = j2 - 1;
          for(int i = 0; i < nx; i++ ) {
              int i2     = g.toIPerm( i );
              int i2_p1  = i2 + 1;
              int i2_n1  = i2 - 1;
              int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
              int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
              int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
              int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
              double cyy_ip1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
              double cyy_in1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
              double cyy_jp1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
              double cyy_jn1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
              int ijk     = g.getIndex( i  , j  , k );
              int ijk_in1 = g.getIndex( i-1, j  , k );
              int ijk_jn1 = g.getIndex( i  , j-1, k );
              permPtr->Perm_x( i, j, k, permx_interface );
              permPtr->Perm_y( i, j, k, permy_interface );
              if(        i == 0      && j == 0) {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );	 
              } else if (i == nx - 1 && j == 0) {
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
              } else if (i ==      0 && j == ny - 1) {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (i == nx - 1 && j == ny - 1) {
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
		 
              } else if (i == 0 ) {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (j == 0 ) {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
              } else if (i == nx - 1) {
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (j == ny - 1) {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else {
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              }
         }
     }
  }
}

void PVar_Sensitivity::calcConstrain(int i1_perm, int j1_perm, 
  int k1_perm, double *solu) {

  double *am1 = new double [nNode];
  double *am2 = new double [nNode];
  double *am3 = new double [nNode];
  double *am4 = new double [nNode];
  for(int i = 0; i < nNode; ++i) am1[i] = am2[i] = am3[i] = am4[i] = 0.;

  double permx_interface[3];
  double permy_interface[3];
  for(int k = 0; k < nz; k++ ) {
      int k2 = g.toKPerm( k );
      for(int j = 0; j < ny; j++ ) {
          int j2 = g.toJPerm( j );
          int j2_p1 = j2 + 1;
          int j2_n1 = j2 - 1;
          for(int i = 0; i < nx; i++ ) {
              int i2     = g.toIPerm( i );
              int i2_p1  = i2 + 1;
              int i2_n1  = i2 - 1;
              int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
              int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
              int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
              int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
              double cyy_ip1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
              double cyy_in1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
              double cyy_jp1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
              double cyy_jn1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
              int ijk     = g.getIndex( i  , j  , k );
              int ijk_in1 = g.getIndex( i-1, j  , k );
              int ijk_jn1 = g.getIndex( i  , j-1, k );
              permPtr->Perm_x( i, j, k, permx_interface );
              permPtr->Perm_y( i, j, k, permy_interface );
              if(        i == 0      && j == 0) {
                 RHS[ijk] = offDAXp[ijk     ] * solu[ijk + 1 ]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );	 
              } else if (i == nx - 1 && j == 0) {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
              } else if (i ==      0 && j == ny - 1) {
                 RHS[ijk] = offDAXp[ijk     ] * solu[ijk + 1 ]
                          + offDAYn[ijk - nx] * solu[ijk - nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (i == nx - 1 && j == ny - 1) {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAYn[ijk - nx] * solu[ijk - nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
		 
              } else if (i == 0 ) {
                 RHS[ijk] = offDAYn[ijk - nx] * solu[ijk - nx]
                          + offDAXp[ijk     ] * solu[ijk + 1 ]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (j == 0 ) {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAXp[ijk     ] * solu[ijk + 1 ]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
              } else if (i == nx - 1) {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAYn[ijk - nx] * solu[ijk - nx]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else if (j == ny - 1) {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAYn[ijk - nx] * solu[ijk - nx]
                          + offDAXp[ijk     ] * solu[ijk + 1 ]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              } else {
                 RHS[ijk] = offDAXn[ijk - 1 ] * solu[ijk - 1 ]
                          + offDAYn[ijk - nx] * solu[ijk - nx]
                          + offDAXp[ijk     ] * solu[ijk + 1 ]
                          + offDAYp[ijk     ] * solu[ijk + nx]
                          +   diagA[ijk     ] * solu[ijk     ];
                 am2[ijk] =-permx_interface[2] * cyy_ip1 
                          * dP0dXi[0][ijk    ] / g.getBdx( i );
                 am1[ijk] = permx_interface[0] * cyy_in1 
                          * dP0dXi[0][ijk_in1] / g.getBdx( i );
                 am4[ijk] =-permy_interface[2] * cyy_jp1 
                          * dP0dXi[1][ijk    ] / g.getBdy( j );
                 am3[ijk] = permy_interface[0] * cyy_jn1 
                          * dP0dXi[1][ijk_jn1] / g.getBdy( j );
              }
         }
     }
  }
 
  for(int i = 0; i < nNode; ++i) {
      if( RHS[i] < 0) {
          RHS[i] = - RHS[i];   	      
          am1[i] = - am1[i];
          am2[i] = - am2[i];
          am3[i] = - am3[i];
          am4[i] = - am4[i];
      }
  }
  
  double tmp, tmp1;
  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(6);
  for(int n =  0; n < nNode; ++n) 
      cout << " res 1 node = " << n <<' '<< RHS[n] << endl;
  
  double upper_bound = 1.;
  double * aa = new double[4];
  double * xx = new double[4];
  for(int n =  0; n < nNode; ++n) {
      if( (n % 2) == 0 ) {
	      
          cout << n <<' ' << RHS[n] <<"   "
               << am1[n] <<' ' << am2[n] <<' '
               << am3[n] <<' ' << am4[n] << endl;
	  aa[0] = am1[n];
	  aa[1] = am2[n];
	  aa[2] = am3[n];
	  aa[3] = am4[n];
	  xx[0] = xx[1] = xx[2] = xx[3] = 0.;
          simplx_search(upper_bound, RHS[n], aa, xx);
	
          int i, j, k = 0;
          j = n / nx;
          i = n - j * nx;
	  if(n != i + j * nx) {
             cout << n << " = " << i + j * nx << endl;
	     exit(0);
	  }

          int i2 = g.toIPerm( i );
          int j2 = g.toJPerm( j );
          int k2 = g.toKPerm( k );
	  /*
	  for(int i = 0; i < 4; ++i) {
              if(xx[i]>1) xx[i] = 1.0;
	  }*/
          permPtr->setYVar(i2 - 1, j2,     k2 , xx[0]);
          permPtr->setYVar(i2 + 1, j2,     k2 , xx[1]);
          permPtr->setYVar(i2    , j2 - 1, k2 , xx[2]);
          permPtr->setYVar(i2    , j2 + 1, k2 , xx[3]);
      }
  }

  delete [] am1;
  delete [] am2;
  delete [] am3;
  delete [] am4;
  delete [] aa;
  delete [] xx;
}

void PVar_Sensitivity::calcCondPStd(int nWells, const Well *w, 
                const Control &c, Solver &s, int iter ) {
  
  ( (P0Eqn*) P0e )->solve(nWells, w, c, s, 0);
      
  // --- point to P0's derivatives for current timestep -----
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, 0) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, 0) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, 0) ;

  double sum = 0;
  for(int kc = 0; kc < nzPerm; kc++) {
      for(int jc = 0; jc < nyPerm; jc++) {
          for(int ic = 0; ic < nxPerm; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
              std_YY[ijkc] = sqrt( permPtr->getYVar(ic,jc,kc) );
              sum += std_YY[ijkc];
          }
      }
  }
  sum /= (nzPerm * nxPerm * nyPerm);
  if(debug)
     cout << "global conditonal std = " <<sum <<endl;

  int kc = 0;
  int jc = 0;
  int ic = 0;
  int ijkc = g.getPermIndex(ic, jc, kc);
  for(int i = 0; i < nNode; ++i) {
      if(fabs(std_PP_UnCd[i]) > 0.00001) {
         rho_YP[i] = CYPe->getCYP(i,ijkc)
                     / std_PP_UnCd[i]/ std_YY_UnCd[ijkc];
      }
      else {
         rho_YP[i] = 0.0;
      }
  } 
  
  calcARHO(c.getBType(), c.getBCond(), nWells, w, 2);
  calcRHS(ic, jc, kc, c.getBType(), c.getBCond());
  s.solve( true, *this, 1, RHS );
  /*
  for(int i = 0; i < nNode; ++i) {
      cout << offDAXn[i-1 ] <<' '
           << offDAXp[i   ] <<' '
           << offDAYn[i-nx] <<' '
           << offDAYp[i   ] <<' '
           <<   diagA[i   ] <<' '
           <<     RHS[i   ] << endl;
  }
  cout << "exit at Pvar()" << endl;
  exit(0);
  */
  /*
  int itmp  = 189;
  cout << "The original RHS = " << RHS[ itmp ] << endl;
  double tmp_left  = offDAXn[itmp - 1 ];
  double tmp_down  = offDAYn[itmp - nx];
  double tmp_right = offDAXp[itmp     ];
  double tmp_upper = offDAYp[itmp     ];
  double tmp_cent =    diagA[itmp     ];

  double  ptmp2 = 0.130061;
  double  ptmp1 = 0.140423;
  
  double  dptmp = 0.0001;

  RHS[itmp] =  ( tmp_left  + tmp_down ) * ptmp1
              +( tmp_right + tmp_upper) * ptmp2;
  
  cout << "node =" << itmp << " rhs before solve= " << RHS[itmp] << endl;
  s.solve( true, *this, 1, RHS );
  cout << "node =" << itmp << " rhs after  solve= " << RHS[itmp] << endl;
  
  double rhs = tmp_left  * RHS[itmp - 1 ]
               + tmp_down  * RHS[itmp - nx]
               + tmp_right * RHS[itmp + 1 ]
             + tmp_upper * RHS[itmp + nx]
             + tmp_cent  * RHS[itmp     ];
  cout << "check: node =" << itmp << " rhs= " << rhs << endl;
  exit(0);
  cout << "Starting Looping " << endl;
  int count = 0;
  while( fabs( RHS[itmp] ) > 0.001 ) {
      double solu_old = RHS[itmp];
      
      calcARHO(c.getBType(), c.getBCond(), nWells, w, 2);
      calcRHS(ic, jc, kc, c.getBType(), c.getBCond());
      
      ptmp1 += dptmp;
      ptmp2 += dptmp;
      RHS[itmp] =  ( tmp_left  + tmp_down ) * ptmp1
                       +( tmp_right + tmp_upper) * ptmp2;
      
      //cout << "node =" << itmp << " rhs= " << RHS[itmp] << endl;
      s.solve( true, *this, 1, RHS );
      cout << "!! ptmp = " << ptmp1 << ' '<<ptmp2 
           << "Loop: solut  = " << RHS[itmp] << " !!!" << endl;
      count++;
      if(count == 5000) exit(0);
      
      //if( (RHS[itmp] - solu_old) < 0) dptmp = -(dptmp);
      //if( fabs( RHS[itmp] ) > fabs(solu_old) ) dptmp /= 2;
  }
  cout << "!! ptmp = " << ptmp1 << ' '<<ptmp2
       << " solut  = " << RHS[itmp] << " !!!" << endl;

  exit(0); 
  */
  for(int i = 0; i < nNode; ++i) {
      if(fabs(RHS[i]) > 0.0000001)
         std_PP[i] = RHS[i];
      else 
         std_PP[i] = 0.0; 
  }
  cout << "std_PP is updated in Pvar" << endl;
}

//===== solve() =====
void PVar_Sensitivity::solve(int nWells, const Well *w, const Control &c, 
                             Solver &s, int iter ) {
                             
  // --- point to P0's derivatives for current timestep -------
  if( nx > 1 ) dP0dXi[0] = P0e->getDP0DXi( 0, 0) ;
  if( ny > 1 ) dP0dXi[1] = P0e->getDP0DXi( 1, 0) ;
  if( nz > 1 ) dP0dXi[2] = P0e->getDP0DXi( 2, 0) ;
  
  //ofstream os("PVar_CPP.out", ios::out);

  int calcType = 2;
  for(int kc = 0; kc < 1; kc++) {
      for(int jc = 0; jc < 1; jc++) {
          for(int ic = 0; ic < 1; ic++) { 
              int ijkc = g.getPermIndex(ic, jc, kc);
              
              for(int i = 0; i < nNode; ++i) {
                  if(fabs(std_PP[i]) > 0.00001) {
                     rho_YP[i] = CYPe->getCYP(i,ijkc)
                                 / std_PP[i]/ std_YY[ijkc];
                  }
                  else {
                     rho_YP[i] = 0.0;
                  }
              } 
              
              calcARHO(c.getBType(), c.getBCond(), nWells, w, calcType);
              calcRHS(ic, jc, kc, c.getBType(), c.getBCond());
              s.solve( true, *this, 1, RHS );
              /*
              if(calcType == 1) {
                 for(int i = 0; i < nNode; ++i) {
                     os << RHS[i] << ' ' << CYPe->getCYP(i,ijkc)<<endl;
                     if( ((i+1)%nx) ==0 ) os << endl;
                 } 
              }
              else if(calcType == 2) {
                 for(int i = 0; i < nNode; ++i) {
                     os << RHS[i] << ' ' << std_PP[i]<<endl;
                     if( ((i+1)%nx) ==0 ) os << endl;
                     if(fabs(RHS[i] - std_PP[i]) > 0.00001) {
                        cout << "wrong! : "<<ic<<' ' << jc<<' '
                             << kc<<' ' <<i <<' '
                             << RHS[i] << ' ' << std_PP[i]<<endl;
                        exit(0);
                     }
                 } 
              }
              else if(calcType == 3) {
                 for(int i = 0; i < nNode; ++i) {
                     os << RHS[i] << ' ' << rho_YP[i]<<endl;
                     if( ((i+1)%nx) ==0 ) os << endl;
                     
                     if(fabs(RHS[i] - rho_YP[i]) > 0.00001) {
                        cout << "wrong! : "<<ic<<' ' << jc<<' '
                             << kc<<' ' <<i <<' '
                             << RHS[i] << ' ' << std_PP[i]<<endl;
                        exit(0);
                     }
                     
                 } 
              }
              */
          }
      }
  }
 // os.close();
  
  if(debug) output();
  cout << "Done at PVar_Sensitivity::solve () " <<endl;
}

//============ Calculate Matrix A ===============
void PVar_Sensitivity::calcARHO(const BType *bType, const double *bCond,
                int nWells, const Well *w, int calcType) {
  int i, j, k, ijk, ijk_p1, ijk_n1;
  double tmp_p, tmp_n, dd;
  double perm_interface[3];
                
  // --- set A to zero ----------
  for(i = 0; i < nNode; i++) diagA[i] = 0.0;
  if( nx > 1 ) 
      for(i = 0; i < nNode - 1; i++) 
          offDAXp[i] = offDAXn[i] = 0.0;
  if( ny > 1 ) 
      for(i = 0; i < nNode - nx; i++) 
          offDAYp[i] = offDAYn[i] = 0.0;
  if( nz > 1 ) 
      for(i = 0; i < nNode - nx*ny; i++) 
          offDAZp[i] = offDAZn[i] = 0.0;
                   
  // --- X direction ----------------------
  if( nx > 1 ) {
      for(k = 0; k < nz; k++ ) {
          for(j = 0; j < ny; j++ ) {
              for(i = 1; i < nx - 1; i++ ) {
                  ijk_p1 = g.getIndex(i + 1, j, k );
                  ijk    = g.getIndex(i    , j, k );
                  ijk_n1 = g.getIndex(i - 1, j, k );
                  permPtr->Perm_x( i, j, k, perm_interface );
                  if(calcType == 1) {
                     offDAXp[ijk  ] = perm_interface[2]  
                                      / g.getBdx( i ) / g.getDx( i );
                     offDAXn[ijk-1] = perm_interface[0] 
                                      / g.getBdx( i ) / g.getDx( i - 1);
                       diagA[ijk ] -= (  perm_interface[2] / g.getDx( i )  
                                       + perm_interface[0] / g.getDx( i - 1)
                                      ) / g.getBdx( i ) ;
                  } else if (calcType == 2) {
                     offDAXp[ijk  ] = perm_interface[2] * rho_YP[ijk_p1] 
                                      / g.getBdx( i ) / g.getDx( i );
                     offDAXn[ijk-1] = perm_interface[0] * rho_YP[ijk_n1]
                                      / g.getBdx( i ) / g.getDx( i - 1);
                       diagA[ijk ] -= (  perm_interface[2] / g.getDx( i    )  
                                       + perm_interface[0] / g.getDx( i - 1)
                                      ) * rho_YP[ijk] / g.getBdx( i ) ;
                    //cout<< " ijk    = "<<ijk
                    //    << " ijk_p1 = "<<ijk_p1
                    //    << " ijk_n1 = "<<ijk_n1 << endl;
                  } else if (calcType == 3) {  
                     offDAXp[ijk  ] = perm_interface[2] * std_PP[ijk_p1] 
                                      / g.getBdx( i ) / g.getDx( i );
                     offDAXn[ijk-1] = perm_interface[0] * std_PP[ijk_n1]
                                      / g.getBdx( i ) / g.getDx( i - 1);
                       diagA[ijk ] -= (  perm_interface[2] / g.getDx( i )  
                                       + perm_interface[0] / g.getDx( i - 1)
                                      ) * std_PP[ijk] / g.getBdx( i );
                  } else {
                    cout << "no implement" << endl;
                    exit(0);
                  }
              } 
          }
      }
      if(g.getGridType() == POINT_DIS ) {
         i      = 0 ;
         double gdx  = g.getDx( 0),  gdx2 = gdx  * gdx;
         double gbdx = g.getBdx( 0), gbdx2 = gbdx * gbdx;
         double tmp;
         switch( bType[0] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) {
                     ijk = g.getIndex(i, j, k);
                     if( j < ny - 1 ) offDAYp[ijk      ] = 0.0;
                     if( j > 0)       offDAYn[ijk - nx ] = 0.0;
                     if( k < nz - 1 ) offDAZp[ijk      ] = 0.0;
                     if( k > 0 )      offDAZn[ijk-nx*ny] = 0.0;
                                        diagA[ijk      ] = 1.0;
                 }
             }
             break;
           case CONST_RATE:
             tmp = 2./gdx2;     
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) {
                     ijk = g.getIndex(0, j, k); 
                     permPtr->Perm_x( 0, j, k, perm_interface );
                     offDAXp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAXp[ijk];
                 }
             }
             break;  
         }
         i    = nx - 1 ;
         gdx  = g.getDx( nx - 2 ),  gdx2 = gdx  * gdx;
         gbdx = g.getBdx( nx - 1 ), gbdx2 = gbdx * gbdx;
         switch( bType[1] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) {
                     ijk = g.getIndex(nx - 1, j, k); 
                     if( j < ny - 1) offDAYp[ijk]       = 0.0;
                     if( j > 0)      offDAYn[ijk-nx]    = 0.0;
                     if( k < nz - 1) offDAZp[ijk]       = 0.0;
                     if( k > 0)      offDAZn[ijk-nx*ny] = 0.0;
                                      diagA[ ijk ]     = 1.0;
                 } 
             }
             break;
           case CONST_RATE:
             tmp = 2./gdx2;
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) {
                     ijk = g.getIndex(nx - 1, j, k); 
                     permPtr->Perm_x( nx - 1, j, k, perm_interface );
                     offDAXn[ijk-1] = tmp * perm_interface[0];
                       diagA[ijk]  -= offDAXn[ijk-1];
                 }
             }
             break;  
         }
      } else {
         i      = 0 ;
         double gdx  = g.getDx( 0),  gdx2 = gdx  * gdx;
         double gbdx = g.getBdx( 0), gbdx2 = gbdx * gbdx;
         double tmp = 1.0/gbdx/gdx;
         switch( bType[0] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) {
                     ijk = g.getIndex(0, j, k);
                     permPtr->Perm_x( 0, j, k, perm_interface );
                     offDAXp[ijk]  = tmp * perm_interface[2]; 
                       diagA[ijk] -= offDAXp[ijk] 
                                   + 2. * perm_interface[0]/gbdx2; 
                 }
             }
             break;
           case CONST_RATE:  // tested in this part
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) { 
                     permPtr->Perm_x( i, j, k, perm_interface );
                     ijk_p1 = g.getIndex(i + 1, j, k);
                     ijk    = g.getIndex(i    , j, k);
                     if(calcType == 1) {
                        offDAXp[ijk]   = tmp * perm_interface[2];
                        diagA[ijk]  -= tmp * perm_interface[2];
                     }
                     else if(calcType == 2) {
                        offDAXp[ijk]   = tmp * perm_interface[2] * rho_YP[ijk_p1];
                          diagA[ijk]  -= tmp * perm_interface[2] * rho_YP[ijk   ];
                     }  
                     else if(calcType == 3) {
                        offDAXp[ijk]   = tmp * perm_interface[2] * std_PP[ijk_p1];
                          diagA[ijk]  -= tmp * perm_interface[2] * std_PP[ijk   ];
                     }
                 }
             }
             break; 
         }
         i      = nx - 1;
         gdx  = g.getDx( nx - 2 ),  gdx2 = gdx  * gdx;
         gbdx = g.getBdx( nx - 1 ), gbdx2 = gbdx * gbdx;
         tmp = 1.0/gbdx/gdx;
         switch( bType[1] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) { 
                     ijk = g.getIndex(nx - 1, j, k); 
                     permPtr->Perm_x( nx - 1, j, k, perm_interface );
                     offDAXn[ijk-1] = tmp * perm_interface[0];
                        diagA[ijk] -= offDAXn[ijk-1] 
                                    + 2. * perm_interface[2]/gbdx2;
                 }
             }
             break;
           case CONST_RATE:  // tested in this part
             for(k = 0; k < nz; k++) {
                 for(j = 0; j < ny; j++) { 
                     permPtr->Perm_x( i, j, k, perm_interface );
                     ijk    = g.getIndex(i    , j, k);
                     ijk_n1 = g.getIndex(i - 1, j, k);
                     if(calcType == 1) {
                        offDAXn[ijk-1] = tmp * perm_interface[0];
                          diagA[ijk ] -= tmp * perm_interface[0];
                     } else if(calcType == 2) {
                        offDAXn[ijk-1] = tmp * perm_interface[0] * rho_YP[ijk_n1];
                          diagA[ijk ] -= tmp * perm_interface[0] * rho_YP[ijk   ];
                          /*
                          cout << offDAXn[ijk-1] <<' '
                               << offDAXp[ijk  ] <<' '
                               << offDAYn[ijk-nx]<<' '
                               << offDAYp[ijk]   <<' '
                               <<   diagA[ijk ]  <<' '
                               <<   RHS[ijk]     <<< endl;
                          exit(0);
                          */
                          //cout << "PVAR:rho_YP = " << rho_YP[ijk_n1] << ' '
                          //     << "PVAR:rho_YP = " << rho_YP[ijk   ] << endl;
                          //exit(0);
                     } else if(calcType == 3) {
                        offDAXn[ijk-1] = tmp * perm_interface[0] * std_PP[ijk_n1];
                          diagA[ijk ] -= tmp * perm_interface[0] * std_PP[ijk   ];
                     } 
                 }
             }
             break;
         }
      }
  }
  
  // --- Y direction ----------------------
  if( ny > 1) {              
      for(j = 1; j < ny - 1; j++) {
          for(k = 0; k < nz; k++) {
              for(i = 0; i < nx; i++) {
                  ijk_p1 = g.getIndex(i, j + 1, k );
                  ijk    = g.getIndex(i, j    , k );
                  ijk_n1 = g.getIndex(i, j - 1, k );
                  permPtr->Perm_y( i, j, k, perm_interface );
                  if(calcType == 1) {
                     offDAYp[ijk   ] = perm_interface[2] 
                                    / g.getBdy( j ) / g.getDy( j  );
                     offDAYn[ijk-nx] = perm_interface[0] 
                                    / g.getBdy( j ) / g.getDy( j - 1);
                       diagA[ijk] -= (  perm_interface[2] / g.getDy( j )
                                      +perm_interface[0] / g.getDy( j - 1) 
                                    ) / g.getBdy( j );
                  }
                  else if(calcType == 2) {
                     offDAYp[ijk   ] = perm_interface[2] * rho_YP[ijk_p1]
                                       / g.getBdy( j ) / g.getDy( j  );
                     offDAYn[ijk-nx] = perm_interface[0] * rho_YP[ijk_n1]
                                       / g.getBdy( j ) / g.getDy( j - 1);
                       diagA[ijk   ]-= ( perm_interface[2] / g.getDy( j )
                                       + perm_interface[0] / g.getDy( j - 1) 
                                       ) * rho_YP[ijk] / g.getBdy( j );
                    //cout<< " ijk    = "<< ijk
                    //    << " ijk_p1 = "<< ijk_p1
                    //    << " ijk_n1 = "<< ijk_n1 << endl;
                  }
                  else if(calcType == 3) {
                     offDAYp[ijk   ] = perm_interface[2] * std_PP[ijk_p1]
                                    / g.getBdy( j ) / g.getDy( j  );
                     offDAYn[ijk-nx] = perm_interface[0] * std_PP[ijk_n1]
                                    / g.getBdy( j ) / g.getDy( j - 1);
                       diagA[ijk] -= (  perm_interface[2] / g.getDy( j )
                                      +perm_interface[0] / g.getDy( j - 1) 
                                    ) * std_PP[ijk] / g.getBdy( j );
                  }
              }                               
          }                                       
      }
      if(g.getGridType() == POINT_DIS ) {
         j      = 0 ;
         double gdy  = g.getDy( 0 ),  gdy2 = gdy  * gdy;
         double gbdy = g.getBdy( 0 ), gbdy2 = gbdy * gbdy;
         double tmp;
         switch( bType[2] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, 0, k);
                     if(i < nx-1) offDAXp[ijk      ] = 0.0;
                     if(i > 0   ) offDAXn[ijk-1    ] = 0.0;
                     if(k < nz-1) offDAZp[ijk      ] = 0.0;
                     if(k > 0   ) offDAZn[ijk-nx*ny] = 0.0;
                                   diagA[ijk      ] = 1.0;
                 }
             }
             break;
           case CONST_RATE:
             tmp = 2./gdy2;
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, 0, k);
                     permPtr->Perm_y( i, 0, k, perm_interface );
                     offDAYp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAYp[ijk]; 
                 }
             }
             break;  
         }
         j    = ny - 1 ;
         gdy  = g.getDy( ny - 2 ),  gdy2 = gdy  * gdy;
         gbdy = g.getBdy( ny - 1 ), gbdy2 = gbdy * gbdy;
         switch( bType[3] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, ny - 1, k);
                     if( i < nx - 1) offDAXp[ijk]   = 0.0;
                     if( i > 0     ) offDAXn[ijk-1] = 0.0;
                     if( k < nz - 1) offDAZp[ijk]   = 0.0;
                     if( k > 0     ) offDAZn[ijk-1] = 0.0;
                                      diagA[ijk]   = 1.0;
                 }
             }
             break;
           case CONST_RATE:
             tmp = 2./gdy2;
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, ny - 1, k);
                     permPtr->Perm_y( i, ny - 1, k, perm_interface );      
                     offDAYn[ijk-nx] = tmp * perm_interface[0];
                       diagA[ijk]   -= offDAYn[ijk-nx]; 
                 }
             }
             break;  
         }
      } else {
         j      = 0 ;
         double gdy  = g.getDy( 0 ),  gdy2 = gdy  * gdy;
         double gbdy = g.getBdy( 0 ), gbdy2 = gbdy * gbdy;
         double tmp = 1.0/gbdy/gdy;
         switch( bType[2] ) {
           case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) { 
                     ijk = g.getIndex(i, 0, k);
                     permPtr->Perm_y( i, 0, k, perm_interface ); 
                     offDAYp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAYp[ijk] 
                                   + 2. * perm_interface[0]/gbdy2;
                 }
             }
             break;
           case CONST_RATE:  // tested in this part
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) { 
                     permPtr->Perm_y( i, j, k, perm_interface ); 
                     ijk_p1 = g.getIndex(i, j + 1, k);
                     ijk    = g.getIndex(i, j    , k);
                     if(calcType == 1) {             
                        offDAYp[ijk]  = tmp * perm_interface[2];
                          diagA[ijk] -= tmp * perm_interface[2];
                     }
                     else if(calcType == 2) { 
                        offDAYp[ijk]  = tmp * perm_interface[2] * rho_YP[ijk_p1];
                          diagA[ijk] -= tmp * perm_interface[2] * rho_YP[ijk   ];
                     }
                     else if(calcType == 3) { 
                        offDAYp[ijk]  = tmp * perm_interface[2] * std_PP[ijk_p1];
                          diagA[ijk] -= tmp * perm_interface[2] * std_PP[ijk   ];
                     }
                 }
             }
             break;  
         }
         j    = ny - 1 ;
         gdy  = g.getDy( ny - 2 ),  gdy2 = gdy  * gdy;
         gbdy = g.getBdy( ny - 1 ), gbdy2 = gbdy * gbdy;
         tmp = 1.0/gbdy/gdy;
         switch( bType[3] ) {
             case CONST_PRES:
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) { 
                     ijk = g.getIndex(i, ny - 1, k); 
                     permPtr->Perm_y( i, ny - 1, k, perm_interface );
                     offDAYn[ijk-nx] = tmp * perm_interface[0];
                       diagA[ijk  ] -= offDAYn[ijk-nx] 
                                     + 2. * perm_interface[2]/gbdy2;
                 }
             }
             break;
           case CONST_RATE:  // tested in this part
             for(k = 0; k < nz; k++) {
                 for(i = 0; i < nx; i++) { 
                     permPtr->Perm_y( i, j, k, perm_interface );
                     ijk    = g.getIndex(i, j    , k);
                     ijk_n1 = g.getIndex(i, j - 1, k);
                     if(calcType == 1) { 
                        offDAYn[ijk-nx] = tmp * perm_interface[0];
                          diagA[ijk   ]-= tmp * perm_interface[0];
                     }
                     else if(calcType == 2) { 
                        offDAYn[ijk-nx] = tmp * perm_interface[0] * rho_YP[ijk_n1];
                          diagA[ijk   ]-= tmp * perm_interface[0] * rho_YP[ijk   ];
                          /*
                          cout << offDAYn[ijk-nx]<<' '
                               << offDAYp[ijk]   <<' '
                               <<   diagA[ijk ]    << endl;
                          exit(0);
                          */
                     }
                     else if(calcType == 3) { 
                        offDAYn[ijk-nx] = tmp * perm_interface[0] * std_PP[ijk_n1];
                          diagA[ijk   ]-= tmp * perm_interface[0] * std_PP[ijk   ];
                     }
                 }
             } 
             break;  
         }
      } 
  }                                        
                                                  
  // --- Z direction (include gravity)--------------
  if( nz > 1 ) {                      
      for(k = 1; k < nz - 1; k++) {
          for(j = 0; j < ny; j++) {
              for(i = 0; i < nx; i++) {
                  ijk_p1 = g.getIndex(i, j, k + 1);
                  ijk    = g.getIndex(i, j, k    );
                  ijk_n1 = g.getIndex(i, j, k - 1);
                  ijk = g.getIndex( i, j, k ); 
                  permPtr->Perm_z( i, j, k, perm_interface );
                  offDAZp[ijk      ] = perm_interface[2] * rho_YP[ijk_p1]
                                       / g.getBdz(k) / g.getDz(k);
                  offDAZn[ijk-nx*ny] = perm_interface[0] * rho_YP[ijk_n1]
                                       / g.getBdz(k) / g.getDz(k-1);
                    diagA[ijk]      -= (  perm_interface[2] / g.getDz(k)
                                        + perm_interface[0] / g.getDz(k-1)
                                       ) * rho_YP[ijk] / g.getBdz(k) ;
                  offDAZp[ijk      ] += perm_interface[1]*fluidPtr->getDen_cf(ijk)
                                        / 2.0 /g.getBdz(k);
                  offDAZn[ijk-nx*ny] += -perm_interface[1]*fluidPtr->getDen_cf(ijk)
                                        / 2.0 /g.getBdz(k);
                  RHS[ijk] = -fluidPtr->getGamma(ijk) *
                            ( perm_interface[2] - perm_interface[0] ) 
                            / g.getBdz(k);
              }                     
          }                       
      }
      if(g.getGridType() == POINT_DIS ) {
         k      = 0 ;
         double gdz  = g.getDz( 0 ),  gdz2 = gdz  * gdz;
         double gbdz = g.getBdz( 0 ), gbdz2 = gbdz * gbdz;
         double tmp;
         switch( bType[4] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, j, 0); 
                     if(i < nx-1) offDAXp[ijk   ] = 0.0;
                     if(i > 0   ) offDAXn[ijk-1 ] = 0.0;
                     if(j < ny-1) offDAYp[ijk   ] = 0.0;
                     if(j > 0   ) offDAYn[ijk-nx] = 0.0;
                                    diagA[ijk   ] = 1.0; 
                 }
             }
             break;
           case CONST_RATE:
             tmp = 2./gdz2;     
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, j, 0);
                     permPtr->Perm_z( i, j, 0, perm_interface );
                     offDAZp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAZp[ijk];
                 }
             }
             break;  
         }
         k    = nz - 1 ;
         gdz  = g.getDz( nz - 2 ),  gdz2 = gdz  * gdz;
         gbdz = g.getBdz( nz - 1 ), gbdz2 = gbdz * gbdz;
         switch( bType[5] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, j, nz - 1); 
                     if( i < nx-1) offDAXp[ijk   ] = 0.0;
                     if( i > 0   ) offDAXn[ijk-1 ] = 0.0;
                     if( j < ny-1) offDAYp[ijk   ] = 0.0;
                     if( j > 0   ) offDAYn[ijk-nx] = 0.0;
                                     diagA[ijk   ] = 1.0;
                 }
             }
             break;
           case CONST_RATE:
             tmp = 2./gdz2;
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, j, nz - 1);
                     permPtr->Perm_z( i, j, nz - 1, perm_interface );
                     offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                       diagA[ijk]      -= offDAZn[ijk-nx*ny];
                 }
             }
             break;  
         }
       } else {
         k      = 0 ;
         double gdz  = g.getDz( 0 ),  gdz2 = gdz  * gdz;
         double gbdz = g.getBdz( 0 ), gbdz2 = gbdz * gbdz;
         double tmp = 1.0/gbdz/gdz;
         switch( bType[4] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) { 
                     ijk = g.getIndex(i, j, 0); 
                     permPtr->Perm_z( i, j, 0, perm_interface );
                     offDAZp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAZp[ijk] 
                                   + 2. * perm_interface[0]/gbdz2;
                 }
             }
             break;
           case CONST_RATE:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) { 
                     ijk = g.getIndex(i, j, 0); 
                     permPtr->Perm_z( i, j, 0, perm_interface );
                     offDAZp[ijk]  = tmp * perm_interface[2];
                       diagA[ijk] -= offDAZp[ijk];
                 }
             }
             break;  
         }
         k    = nz - 1 ;
         gdz  = g.getDz( nz - 2 ),  gdz2 = gdz  * gdz;
         gbdz = g.getBdz( nz - 1 ), gbdz2 = gbdz * gbdz;
         tmp  = 1.0/gbdz/gdz;
         switch( bType[5] ) {
           case CONST_PRES:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) {
                     ijk = g.getIndex(i, j, nz - 1);
                     permPtr->Perm_z( i, j, nz - 1, perm_interface );
                     offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                       diagA[ijk]      -= offDAZn[ijk-nx*ny] 
                                         + 2. * perm_interface[2]/gbdz2;
                 }
             }
             break;
           case CONST_RATE:
             for(j = 0; j < ny; j++) {
                 for(i = 0; i < nx; i++) { 
                     ijk = g.getIndex(i, j, nz - 1);
                     permPtr->Perm_z( i, j, nz - 1, perm_interface );
                     offDAZn[ijk-nx*ny] = tmp * perm_interface[0];
                       diagA[ijk]      -= offDAZn[ijk-nx*ny]; 
                 }
             }
             break;  
         }
     } 
  }
  // well
  for(int iw = 0; iw < nWells; iw++ ) {
      if(w[iw].getWellType() == PRES) {
         for(int i = 0; i < w[iw].getWellLen(); i++ ) {
             int ii  = w[iw].getWellIJK(0, i);
             int jj  = w[iw].getWellIJK(1, i);
             int kk  = w[iw].getWellIJK(2, i);      
             int ijk = g.getIndex(ii, jj, kk); 
             double vol  = g.getBdx( ii ) * g.getBdy( jj ) * g.getBdz( kk );
             diagA[ijk] -=  w[iw].getWellTrans(i)/vol;
         }
      }
  }
}

void PVar_Sensitivity::setCondPress() {
  for(int m = 0; m < num_meas; m++ ) {
      int i   = prsmPtr->getIpmeas(m);
      int j   = prsmPtr->getJpmeas(m);
      int k   = prsmPtr->getKpmeas(m);
      int ijk = g.getIndex(i, j, k); 
      double vol = g.getBdx( i ) * g.getBdy( j ) * g.getBdz( k );
      diagA[ijk] -=  bigNumber/vol;
    //  cout <<"i = "
    //           << i << ' '<<j<<' '<<k <<' '<<ijk<<' '<<diagA[ijk]<<endl;
  }
  //exit(0);
}

//===== ( calcSensiRHS() ) =====
void PVar_Sensitivity::calcSensiRHS(int& master, int & i1_perm, int & j1_perm, 
                 int & k1_perm, const BType *bType, const double *bCond) {
                                
  int nxy = nx * ny;
  int i , j , k ;
  int i2, j2, k2;
  int ijk, ijk_n1;
  int i2_p1, i2_n1;
  int j2_p1, j2_n1;
  int k2_p1, k2_n1;

  double   cyy_p1, cyy_n1, cyy_p0;
  double perm_interface[3];
  
  int      ijk2_p1,   ijk2_n1,   ijk2_p0;
  double   lambda_p1, lambda_n1, lambda_p0; 
  
  double std_Y_master = std_YY[g.getPermIndex(i1_perm, j1_perm, k1_perm)];
  
  for( i = 0; i < nNode; i++ ) RHS[i] = 0.0;
  
  if(nx > 1 ) {
     double dx0   = g.getDx(0);
     double dx02  = dx0 * dx0;
     double dxn   = g.getDx( nx - 2 );
     double gbdx0 = g.getBdx( 0 );
     double gbdxn = g.getBdx( nx - 1);
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(j = 0; j < ny; j++ ) {
             j2 = g.toJPerm( j );
             for(i = 1; i < nx - 1; i++ ) {
                 i2     = g.toIPerm( i );

                 i2_p1  = i2 + 1;
                 i2_n1  = i2 - 1;
                 ijk2_p1 = g.getPermIndex(i2_p1, j2, k2);
                 ijk2_n1 = g.getPermIndex(i2_n1, j2, k2);
                 lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                 lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                 
                 ijk_n1 = g.getIndex( i-1, j, k );
                 ijk    = g.getIndex( i  , j, k );
                 permPtr->Perm_x( i, j, k, perm_interface );
                 RHS[ijk] -= ( perm_interface[2] * lambda_p1 * cyy_p1 * dP0dXi[0][ijk] -
                               perm_interface[0] * lambda_n1 * cyy_n1 * dP0dXi[0][ijk_n1]
                             ) / g.getBdx( i ); 
                 /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<<"nx>1"<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/
             }
         }
     }

     if(g.getGridType() == POINT_DIS ) {
        i     = 0 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2 = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*( cyy_p1 * dP0dXi[0][ijk]
                                   + bCond[0] * (cyy_p1 - cyy_p0)/perm_interface[1]
                                   ) / dx0 * perm_interface[2];
                }
            }
            break;
        }
        i     = nx - 1 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++) 
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++) {
                    j2     = g.toJPerm( j );
                    ijk    = g.getIndex(     i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dxn 
                                  * ( cyy_n1 * ( -dP0dXi[0][ijk_n1] )
                                    - bCond[1] * (cyy_p0 - cyy_n1 )/ perm_interface[1]
                                    );
                }
            }
            break;
        }
      }
     else {
        i     = 0 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[0])/(0.5*gbdx0)
                                ) / gbdx0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );

                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    ijk2_p1    = g.getPermIndex(i2_p1, j2, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                    
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * lambda_p1 * cyy_p1 
                               / gbdx0 * dP0dXi[0][ijk];
                    /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/

                }
            }
            break;
        }
        i     = nx - 1;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:  
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                   (bCond[1] - p0[ijk])/(0.5 * gbdxn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[0][ijk_n1])
                                ) / gbdxn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );

                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk2_n1    = g.getPermIndex(i2_n1, j2, k2);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
        
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyy_n1 * lambda_n1
                                / gbdxn * (-dP0dXi[0][ijk_n1]);
                    /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/

                }
            }
            break;
        }
     }
  }

  if(ny > 1 ) {
     double dy0   = g.getDy( 0 );
     double dy02  = dy0 * dy0;
     double dyn   = g.getDy( ny - 2 );
     double gbdy0 = g.getBdy( 0 );
     double gbdyn = g.getBdy( ny - 1 );
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(j = 1; j < ny - 1; j++ ) {
                 j2     = g.toJPerm( j );
                 j2_p1  = j2 + 1;
                 j2_n1  = j2 - 1;
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);

                 ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                 ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                 lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);
                 lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);
                 
                 ijk_n1 = g.getIndex( i, j-1, k );
                 ijk    = g.getIndex( i, j  , k );
                 permPtr->Perm_y( i, j, k, perm_interface );
                 RHS[ijk] -= (  perm_interface[2] * cyy_p1 * lambda_p1 * dP0dXi[1][ijk ]
                              - perm_interface[0] * cyy_n1 * lambda_n1 * dP0dXi[1][ijk_n1 ]
                             ) / g.getBdy( j );
                 /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/
             }
         }
     }

     if(g.getGridType() == POINT_DIS ) {
        j     = 0 ;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        switch( bType[2] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j  , k );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2] /dy0 
                                   * ( cyy_p1 * dP0dXi[1][ijk]
                                      +bCond[2]* (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        j     = ny - 1 ;
        j2    = g.toJPerm( j );
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[0] / dyn
                                   * ( cyy_n1 * ( -dP0dXi[1][ijk_n1] )
                                     - bCond[3]* (cyy_p0 - cyy_n1)/perm_interface[1]
                                     );
                }
            }
            break;
        }
     } 
     else {
        j     = 0 ;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[2] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[2])/(0.5*gbdy0)
                                ) / gbdy0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );

                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    ijk2_p1    = g.getPermIndex(i2, j2_p1, k2);
                    lambda_p1  = krigPtr->getWeigh(master,ijk2_p1);

                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyy_p1 * lambda_p1
                                / gbdy0 * dP0dXi[1][ijk];
                 /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/

                }
            }
            break;
        }
        j     = ny - 1;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES: 
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                    (bCond[3] - p0[ijk])/(0.5*gbdyn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[1][ijk_n1] )
                                ) / gbdyn;
                }
            }
            break;
          case CONST_RATE:  //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );

                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk2_n1    = g.getPermIndex(i2, j2_n1, k2);
                    lambda_n1  = krigPtr->getWeigh(master,ijk2_n1);

                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyy_n1 * lambda_n1
                                / gbdyn * (-dP0dXi[1][ijk_n1]);
                 /*
                 if(fabs(RHS[ijk]) > 100000) {
                   cout<< i <<' ' <<j <<' '<< k<<' ';
                   cout<<  lambda_p1 <<' '<<lambda_n1<<endl;
                   exit(0);
                 }*/

                }
            }
            break;
         }
     }
  }

  if(nz > 1) {
     double dz0   = g.getDz(0);
     double dz02  = dz0 * dz0;
     double dzn   = g.getDz(nz - 2);
     double gbdz0 = g.getBdz(0);
     double gbdzn = g.getBdz(nz - 1);
     for(j = 0; j < ny; j++ ) {
         j2 = g.toJPerm( j );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(k = 1; k < nz - 1; k++ ) {
                 k2 = g.toKPerm( k );
                 k2_p1  = k2 + 1;
                 k2_n1  = k2 - 1;
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                 ijk_n1 = g.getIndex( i, j, k - 1 );
                 ijk    = g.getIndex( i, j, k     );
                 permPtr->Perm_z( i, j, k, perm_interface );
                 RHS[ijk] -=   ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk   ] + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk) )
                               ) / g.getBdz( k );
             }
         }
     }
     if(g.getGridType() == POINT_DIS) {
        k     = 0 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2 = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2]/ dz0 
                                   * ( cyy_p1 * (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                      +bCond[4] * (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        k     = nz - 1 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dzn   
                                  * ( cyy_n1 * ( - dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                    - bCond[5]* (cyy_p0 - cyy_n1 )/perm_interface[1]
                                    );
                }
            }
            break;
        }
     }
     else {
        k     = 0 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:  //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                 -perm_interface[0] * cyy_n1 *
                                  ( (p0[ijk] - bCond[4])/(0.5*gbdz0) + fluidPtr->getGamma(ijk) )
                                ) / gbdz0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                                ) /gbdz0;
                }
            }
            break;
        }
        k     = nz - 1 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                  ((bCond[5] - p0[ijk])/(0.5*gbdzn) + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk))
                                 ) / gbdzn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[0] * cyy_n1 * 
                                  (-dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                 -perm_interface[2] * cyy_p1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdzn;
                }
            }
            break;
        }
     }
  }
}     

//===== ( calcRHS() ) =====
void PVar_Sensitivity::calcRHS( int & i1_perm, int & j1_perm, int & k1_perm, 
                                const BType *bType, const double *bCond) {
  int nxy = nx * ny;
  int i , j , k ;
  int i2, j2, k2;
  int ijk, ijk_p1, ijk_n1;
  int i2_p1, i2_n1;
  int j2_p1, j2_n1;
  int k2_p1, k2_n1;

  double   cyy_p1, cyy_n1, cyy_p0;
  double perm_interface[3];
  
  for( i = 0; i < nNode; i++ ) RHS[i] = 0.0;
  
  // --- X direction (dCYPdx) -------------------------
  if(nx > 1 ) {
     double dx0   = g.getDx(0);
     double dx02  = dx0 * dx0;
     double dxn   = g.getDx( nx - 2 );
     double gbdx0 = g.getBdx( 0 );
     double gbdxn = g.getBdx( nx - 1);
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(j = 0; j < ny; j++ ) {
             j2 = g.toJPerm( j );
             for(i = 1; i < nx - 1; i++ ) {
                 i2     = g.toIPerm( i );
                 i2_p1  = i2 + 1;
                 i2_n1  = i2 - 1;
                 ijk_n1 = g.getIndex( i-1, j, k );
                 ijk    = g.getIndex( i  , j, k );
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                 permPtr->Perm_x( i, j, k, perm_interface );
                 RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk] -
                               perm_interface[0] * cyy_n1 * dP0dXi[0][ijk_n1]
                             ) / g.getBdx( i ); 
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        i     = 0 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2 = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.*( cyy_p1 * dP0dXi[0][ijk]
                                   + bCond[0] * (cyy_p1 - cyy_p0)/perm_interface[1]
                                   ) / dx0 * perm_interface[2];
                }
            }
            break;
        }
        i     = nx - 1 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++) 
                for(j = 0; j < ny; j++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++) {
                    j2     = g.toJPerm( j );
                    ijk    = g.getIndex(     i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2   , j2, k2);
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dxn 
                                  * ( cyy_n1 * ( -dP0dXi[0][ijk_n1] )
                                    - bCond[1] * (cyy_p0 - cyy_n1 )/ perm_interface[1]
                                    );
                }
            }
            break;
        }
      }
     else {
        i     = 0 ;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[0] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[0][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[0])/(0.5*gbdx0)
                                ) / gbdx0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyy_p1 / gbdx0 * dP0dXi[0][ijk];
                }
            }
            break;
        }
        i     = nx - 1;
        i2    = g.toIPerm( i );
        i2_p1 = i2 + 1;
        i2_n1 = i2 - 1;
        switch( bType[1] ) {
          case CONST_PRES:  //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                   (bCond[1] - p0[ijk])/(0.5 * gbdxn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[0][ijk_n1])
                                ) / gbdxn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(j = 0; j < ny; j++ ) {
                    j2     = g.toJPerm( j );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i - 1, j, k );
                    permPtr->Perm_x( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyy_n1 / gbdxn * (-dP0dXi[0][ijk_n1]);
                }
            }
            break;
        }
     }
  }

  if(ny > 1 ) {
     double dy0   = g.getDy( 0 );
     double dy02  = dy0 * dy0;
     double dyn   = g.getDy( ny - 2 );
     double gbdy0 = g.getBdy( 0 );
     double gbdyn = g.getBdy( ny - 1 );
     for(k = 0; k < nz; k++ ) {
         k2 = g.toKPerm( k );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(j = 1; j < ny - 1; j++ ) {
                 j2     = g.toJPerm( j );
                 j2_p1  = j2 + 1;
                 j2_n1  = j2 - 1;
                 ijk_n1 = g.getIndex( i, j-1, k );
                 ijk    = g.getIndex( i, j  , k );
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                 permPtr->Perm_y( i, j, k, perm_interface );
                 RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk ]
                              - 
                               perm_interface[0] * cyy_n1 * dP0dXi[1][ijk_n1 ]
                             ) / g.getBdy( j );
             }
         }
     }
     if(g.getGridType() == POINT_DIS ) {
        j     = 0 ;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        switch( bType[2] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for(i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j  , k );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2] /dy0 
                                   * ( cyy_p1 * dP0dXi[1][ijk]
                                      +bCond[2]* (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        j     = ny - 1 ;
        j2    = g.toJPerm( j );
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES:
            for(k = 0; k < nz; k++)
                for( i = 0; i < nx; i++) RHS[ g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2   , k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[0] / dyn
                                   * ( cyy_n1 * ( -dP0dXi[1][ijk_n1] )
                                     - bCond[3]* (cyy_p0 - cyy_n1)/perm_interface[1]
                                     );
                }
            }
            break;
        }
     } 
     else {
        j     = 0 ;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[2] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * dP0dXi[1][ijk]
                                 -perm_interface[0] * cyy_n1 *
                                    (p0[ijk] - bCond[2])/(0.5*gbdy0)
                                ) / gbdy0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[2] * cyy_p1 / gbdy0 * dP0dXi[1][ijk];
                }
            }
            break;
        }
        j     = ny - 1;
        j2    = g.toJPerm( j );
        j2_p1 = j2 + 1;
        j2_n1 = j2 - 1;
        switch( bType[3] ) {
          case CONST_PRES: //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    ijk    = g.getIndex( i, j    , k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                    (bCond[3] - p0[ijk])/(0.5*gbdyn)
                                 -perm_interface[0] * cyy_n1 * ( dP0dXi[1][ijk_n1] )
                                ) / gbdyn;
                }
            }
            break;
          case CONST_RATE:  //ok
            for(k = 0; k < nz; k++) {
                k2 = g.toKPerm( k );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j - 1, k );
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
                    permPtr->Perm_y( i, j, k, perm_interface );
                    RHS[ijk] -= perm_interface[0] * cyy_n1 / gbdyn * (-dP0dXi[1][ijk_n1]);
                }
            }
            break;
         }
     }
  }

  if(nz > 1) {
     double dz0   = g.getDz(0);
     double dz02  = dz0 * dz0;
     double dzn   = g.getDz(nz - 2);
     double gbdz0 = g.getBdz(0);
     double gbdzn = g.getBdz(nz - 1);
     for(j = 0; j < ny; j++ ) {
         j2 = g.toJPerm( j );
         for(i = 0; i < nx; i++ ) {
             i2 = g.toIPerm( i );
             for(k = 1; k < nz - 1; k++ ) {
                 k2 = g.toKPerm( k );
                 k2_p1  = k2 + 1;
                 k2_n1  = k2 - 1;
                 cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                 cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                 ijk_n1 = g.getIndex( i, j, k - 1 );
                 ijk    = g.getIndex( i, j, k     );
                 permPtr->Perm_z( i, j, k, perm_interface );
                 RHS[ijk] -=   ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk   ] + fluidPtr->getGamma(ijk) )
                                -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk) )
                               ) / g.getBdz( k );
             }
         }
     }
     if(g.getGridType() == POINT_DIS) {
        k     = 0 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2 = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2. * perm_interface[2]/ dz0 
                                   * ( cyy_p1 * (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                      +bCond[4] * (cyy_p1 - cyy_p0 )/perm_interface[1]
                                     );
                }
            }
            break;
        }
        k     = nz - 1 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES:
            for(j = 0; j < ny; j++) 
                for(i = 0; i < nx; i++) RHS[g.getIndex( i, j, k ) ] = 0.0;
            break;
          case CONST_RATE:
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    cyy_p0 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= 2.* perm_interface[0] / dzn   
                                  * ( cyy_n1 * ( - dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                    - bCond[5]* (cyy_p0 - cyy_n1 )/perm_interface[1]
                                    );
                }
            }
            break;
        }
     }
     else {
        k     = 0 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1 = k2 - 1;
        switch( bType[4] ) {
          case CONST_PRES:  //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  (dP0dXi[2][ijk] + fluidPtr->getGamma(ijk) )
                                 -perm_interface[0] * cyy_n1 *
                                  ( (p0[ijk] - bCond[4])/(0.5*gbdz0) + fluidPtr->getGamma(ijk) )
                                ) / gbdz0;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    ijk    = g.getIndex( i,j,k   );
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 * 
                                  ( dP0dXi[2][ijk] + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                                ) /gbdz0;
                }
            }
            break;
        }
        k     = nz - 1 ;
        k2    = g.toKPerm( k );
        k2_p1 = k2 + 1;
        k2_n1  = k2 - 1;
        switch( bType[5] ) {
          case CONST_PRES: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[2] * cyy_p1 *
                                  ((bCond[5] - p0[ijk])/(0.5*gbdzn) + fluidPtr->getGamma(ijk))
                                 -perm_interface[0] * cyy_n1 * 
                                  ( dP0dXi[2][ijk_n1] + fluidPtr->getGamma(ijk))
                                 ) / gbdzn;
                }
            }
            break;
          case CONST_RATE: //ok
            for(j = 0; j < ny; j++) {
                j2 = g.toJPerm( j );
                for(i = 0; i < nx; i++) {
                    i2     = g.toIPerm( i );
                    ijk    = g.getIndex( i, j, k );
                    ijk_n1 = g.getIndex( i, j, k - 1);
                    cyy_p1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_p1);
                    cyy_n1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2, k2_n1);
                    permPtr->Perm_z( i, j, k, perm_interface );
                    RHS[ijk] -= ( perm_interface[0] * cyy_n1 * 
                                  (-dP0dXi[2][ijk_n1] - fluidPtr->getGamma(ijk))
                                 -perm_interface[2] * cyy_p1 *
                                  ( fluidPtr->getGamma(ijk) - fluidPtr->getGamma(ijk) )
                               ) /gbdzn;
                }
            }
            break;
        }
     }
  }
}     

void PVar_Sensitivity::calcAccu(bool recalA, double dt) {
  cout<<"no implementation!" << endl;
}


void PVar_Sensitivity::display() {
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      cout << "masterpoint = " << m << " Sensititity = ";
      for(int i = 0; i < num_meas; ++i) {
          cout << pVar_Sensitivity[ i + m * num_meas] <<' ';
      }
      cout << endl;
  } 
}

void PVar_Sensitivity::output() {
  ofstream os("PVar_sensitivity.out", ios::out);
//ofstream os("PVar_sensitivity.out", ios::app);
  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(4);
  for(int m = 0; m < mptsPtr->getLength(); ++m) {
      os << "m = " << m << ' ';
      for(int i = 0; i < num_meas; ++i) {
          os << pVar_Sensitivity[ i + m * num_meas] <<' ';
      }
      os <<endl;
  } 
  os.close();
}

void PVar_Sensitivity::calcPermStdBySimplx(int & i1_perm, int & j1_perm, int & k1_perm) {
  int icase;
  int n  = nNodePerm;
  n = 4;
  int np = n + 1;   
  int m1 = n;
  int m2 = 0; 
  int m3 = nNode;
  m3 = 3;
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];

  // Initialize to 0.
  for(int i = 0; i < mp * np; ++i) am[i] = 0.0;
  
  // The first row in the tableau, or row 0
  am[ 0 + 0 * mp] =  0.0;
  for(int j = 1; j < np; ++j)  am[ 0 + j * mp] =  1.0;
 
  // Rows from 1 to m1 in the tableau, the part with <=  
  for(int i = 1; i <= m1; ++i) {
      am[ i + 0 * mp] =  60.0;
      am[ i + i * mp] = -1.0;
  }
  
  // 
  int * m3i = new int[m3];
  int * m3j = new int[m3];
  int * m3k = new int[m3];
  m3i[0] = 2;
  m3j[0] = 2;
  m3k[0] = 0;
  int m1p1 = m1 + 1;
  double permx_interface[3];
  double permy_interface[3];
  for(int i = 0; i < m3; ++i) {
      int ijk     = g.getIndex(m3i[i]  , m3j[i]  , m3k[i]);
      int ijk_in1 = g.getIndex(m3i[i]-1, m3j[i]  , m3k[i]);
      int ijk_jn1 = g.getIndex(m3i[i]  , m3j[i]-1, m3k[i]);
      permPtr->Perm_x( m3i[i]  , m3j[i]  , m3k[i], permx_interface );
      permPtr->Perm_y( m3i[i]  , m3j[i]  , m3k[i], permy_interface );

      int k2 = g.toKPerm( m3k[i] );
      int j2 = g.toJPerm( m3j[i] );
      int i2 = g.toIPerm( m3i[i] );
      int j2_p1 = j2 + 1;
      int j2_n1 = j2 - 1;
      int i2_p1 = i2 + 1;
      int i2_n1 = i2 - 1;
      int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
      int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
      int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
      int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
      double cyy_ip1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
      double cyy_in1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
      double cyy_jp1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
      double cyy_jn1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);
      
      am[ i + m1p1 + 0 * mp] = RHS[ijk]; 
      am[ i + m1p1 + 1 * mp] = permx_interface[2] * cyy_ip1 
                             * dP0dXi[0][ijk    ] / g.getBdx( m3i[i] );
      am[ i + m1p1 + 2 * mp] =-permx_interface[0] * cyy_in1 
                             * dP0dXi[0][ijk_in1] / g.getBdx( m3i[i] );
      am[ i + m1p1 + 3 * mp] = permy_interface[2] * cyy_jp1 
                             * dP0dXi[1][ijk    ] / g.getBdy( m3j[i] );
      am[ i + m1p1 + 4 * mp] =-permy_interface[0] * cyy_jn1 
                             * dP0dXi[1][ijk_jn1] / g.getBdy( m3j[i] );
  }
  
  // change sign
  for(int i = 0; i < m3; ++i) {
      if(am[ i + m1p1 + 0 * mp] < 0 ) {
         am[ i + m1p1 + 0 * mp] = - am[ i + m1p1 + 0 * mp];
         for(int j = 1; j < np; ++j) {
             am[ i + m1p1 + j * mp] = -am[ i + m1p1 + j * mp];
         }
      }
  }
 
  for(int i = 0; i < mp; ++i) {
      for(int j = 0; j < np; ++j) {
          cout << am[i + j * mp] << ' ';
      }
      cout << endl;
  }
  //exit(0); 

  /*
  double permx_interface[3];
  double permy_interface[3];
  for(int k = 0; k < nz; k++ ) {
      int k2 = g.toKPerm( k );
      for(int j = 0; j < ny; j++ ) {
          int j2 = g.toJPerm( j );
          int j2_p1 = j2 + 1;
          int j2_n1 = j2 - 1;
          for(int i = 0; i < nx; i++ ) {
              int i2     = g.toIPerm( i );
              int i2_p1  = i2 + 1;
              int i2_n1  = i2 - 1;
              
              int ijk2_ip1 = g.getPermIndex(i2_p1, j2   , k2);
              int ijk2_in1 = g.getPermIndex(i2_n1, j2   , k2);
              int ijk2_jp1 = g.getPermIndex(i2   , j2_p1, k2);
              int ijk2_jn1 = g.getPermIndex(i2   , j2_n1, k2);
              double cyy_ip1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_p1, j2, k2);
              double cyy_in1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2_n1, j2, k2);
              double cyy_jp1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_p1, k2);
              double cyy_jn1 = permPtr->getCYY(i1_perm, j1_perm, k1_perm, i2, j2_n1, k2);

              int ijk     = g.getIndex( i  , j  , k );
              int ijk_in1 = g.getIndex( i-1, j  , k );
              int ijk_jn1 = g.getIndex( i  , j-1, k );
              permPtr->Perm_x( i, j, k, permx_interface );
              permPtr->Perm_y( i, j, k, permy_interface );
              if( ijk < m3 ) {
              if(        i == 0      && j == 0     ) {
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;
              }
              
              else if (i == nx - 1 && j == 0     ) {
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1
                                                 * dP0dXi[0][ijk_in1] / g.getBdx(nx-1);
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (i ==      0 && j == ny - 1) {
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy(ny-1);
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (i == nx - 1 && j == ny - 1) {
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1
                                                 * dP0dXi[0][ijk_in1] / g.getBdx(nx-1);
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy(ny-1); 
                
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (i == 0     ) {
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy( j );
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (j == 0     ) {
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1 
                                                 * dP0dXi[0][ijk_in1] / g.getBdx( i );
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (i == nx - 1) {
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1
                                                 * dP0dXi[0][ijk_in1] / g.getBdx(nx-1);
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy( j );
                
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else if (j == ny - 1) {
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy(ny-1);
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1 
                                                 * dP0dXi[0][ijk_in1] / g.getBdx( i );
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              } 
              
              else {
                am[ijk+np + (ijk2_ip1 + 1) * mp] = permx_interface[2] * cyy_ip1 
                                                 * dP0dXi[0][ijk    ] / g.getBdx( i );
                am[ijk+np + (ijk2_in1 + 1) * mp] =-permx_interface[0] * cyy_in1 
                                                 * dP0dXi[0][ijk_in1] / g.getBdx( i );
                am[ijk+np + (ijk2_jp1 + 1) * mp] = permy_interface[2] * cyy_jp1 
                                                 * dP0dXi[1][ijk    ] / g.getBdy( j );
                am[ijk+np + (ijk2_jn1 + 1) * mp] =-permy_interface[0] * cyy_jn1 
                                                 * dP0dXi[1][ijk_jn1] / g.getBdy( j );
                if( am[ijk+np + 0 * mp] < 0 ) {
                    am[ijk+np + 0 * mp] = -am[ijk+np + 0 * mp];
                    am[ijk+np + (ijk2_ip1 + 1) * mp] = -am[ijk+np + (ijk2_ip1 + 1) * mp];
                    am[ijk+np + (ijk2_in1 + 1) * mp] = -am[ijk+np + (ijk2_in1 + 1) * mp];
                    am[ijk+np + (ijk2_jp1 + 1) * mp] = -am[ijk+np + (ijk2_jp1 + 1) * mp];
                    am[ijk+np + (ijk2_jn1 + 1) * mp] = -am[ijk+np + (ijk2_jn1 + 1) * mp]; 
                }
                cout << ijk << " : "
                     << am[ijk+np + 0 * mp]<< " = "
                     << am[ijk+np + (ijk2_ip1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_in1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jp1 + 1) * mp] << ' '
                     << am[ijk+np + (ijk2_jn1 + 1) * mp] << endl;                
              }
              }
          } 
      }
  }
  */ 
  //exit(0);

  int nmax = n;
  int mmax = m;
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];
  simplx(am, &m, &n, &mp, &np,
         &m1, &m2, &m3, &icase, izrov, iposv,
         &nmax, l1, &mmax, l2, l3);

  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
  } else {
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) 
         cout << "x[" << iposv[i] << "] = "<< am[ i  + 1] << endl;
     }
     cout << endl;
  }

  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

void PVar_Sensitivity::simplx_driver2() {
  int icase;
  int n  = 4;
  int m1 = 2;
  int m2 = 1; 
  int m3 = 1;
  int np = n + 1;  
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];

  int nmax = n;
  int mmax = m;
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];

  am[     0      ] =  0.0;
  am[     1 * mp ] =  1.0;
  am[     2 * mp ] =  1.0;
  am[     3 * mp ] =  3.0;
  am[     4 * mp ] = -0.5;
  
  am[ 1 + 0 * mp ] =740.0;
  am[ 1 + 1 * mp ] = -1.0;
  am[ 1 + 2 * mp ] =  0.0;
  am[ 1 + 3 * mp ] = -2.0;
  am[ 1 + 4 * mp ] =  0.0;
  
  am[ 2 + 0 * mp ] =  0.0;
  am[ 2 + 1 * mp ] =  0.0;
  am[ 2 + 2 * mp ] = -2.0;
  am[ 2 + 3 * mp ] =  0.0;
  am[ 2 + 4 * mp ] =  7.0;
  
  am[ 3 + 0 * mp ] =  0.5;
  am[ 3 + 1 * mp ] =  0.0;
  am[ 3 + 2 * mp ] = -1.0;
  am[ 3 + 3 * mp ] =  1.0;
  am[ 3 + 4 * mp ] = -2.0;

  am[ 4 + 0 * mp ] =  9.0;
  am[ 4 + 1 * mp ] = -1.0;
  am[ 4 + 2 * mp ] = -1.0;
  am[ 4 + 3 * mp ] = -1.0;
  am[ 4 + 4 * mp ] = -1.0;

  for(int j = 0; j < np; ++j) 
      am[ mp - 1 + j * mp] = 0.;
  
  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(4);
  cout << endl; 
  for(int i = 0; i < mp; ++i) {
      for(int j = 0; j < np; ++j) {      
          cout<< am[ i + j * mp] <<' ';
      }
      cout << endl;
  }
  cout << endl;
  
  simplx(am, &m, &n, &mp, &np,
         &m1, &m2, &m3, &icase, izrov, iposv,
         &nmax, l1, &mmax, l2, l3);
  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
  } else {
     for(int i = 0; i < mp; ++i) {
         for(int j = 0; j < np; ++j) {      
             cout<< am[ i + j * mp] <<' ';
         }
         cout << endl;
     }
     cout << endl;
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) 
         cout << "x[" << iposv[i] << "] = "<< am[ i  + 1] << endl;
     }
     cout << endl;
     for(int j = 0; j < n; ++j) cout << izrov[j]<<' ';
     cout<< endl;
  }

  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

void PVar_Sensitivity::simplx_driver1() {
  int icase;
  int n  = 4;
  int m1 = 0;
  int m2 = 0; 
  int m3 = 2;
  int np = n + 1;
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];
  int nmax = n;
  int mmax = m;
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];
  for(int i = 0; i < mp * np; ++i) am[i] = 0;  
  am[1]  = am[8] =2.0;
  am[2]  = 8.0;
  am[10] = 3.0;
  am[13] = 1.0;
  am[9]  = -6.0;
  am[12] = am[14] = -4.0;  
  am[5]  = am[18] = -1.0; 

  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(4);
  cout << endl; 
  for(int i = 0; i < mp; ++i) {
      for(int j = 0; j < np; ++j) {      
          cout<< am[ i + j * mp] <<' ';
      }
      cout << endl;
  }
  cout << endl;
  
  simplx(am, &m, &n, &mp, &np,
         &m1, &m2, &m3, &icase, izrov, iposv,
         &nmax, l1, &mmax, l2, l3);
  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
  } else {
     for(int i = 0; i < mp; ++i) {
         for(int j = 0; j < np; ++j) {      
             cout<< am[ i + j * mp] <<' ';
         }
         cout << endl;
     }
     cout << endl;
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) 
         cout << "x[" << iposv[i] << "] = "<< am[ i  + 1] << endl;
     }
     cout << endl;
     for(int j = 0; j < n; ++j) cout << izrov[j]<<' ';
     cout<< endl;
  }

  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

void PVar_Sensitivity::simplx_search(double upper_bound, double res,
		double* aa, double* x) {
  int icase;
  int n  = 4;
  int m1 = 4;
  int m2 = 0; 
  int m3 = 1;
  int np = n + 1;
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  int nmax = n;
  int mmax = m;
  
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];
  // Initialize to 0.
  for(int i = 0; i < mp * np; ++i) am[i] = 0.0;
  
  // The first row in the tableau, or row 0
  am[ 0 + 0 * mp] =  0.0;
  for(int j = 1; j < np; ++j)  am[ 0 + j * mp] =  1.0;

  // Rows from 1 to m1 in the tableau, the part with <=  
  for(int i = 1; i <= m1; ++i) {
      am[ i + 0 * mp] =  upper_bound;
      am[ i + i * mp] = -1.0;
  }
  
  // matrix 
  am[m1 + 1 + 0 * mp] = res;
  for(int i = 1; i <= n; ++i) {
      am[m1 + 1 + i * mp] = -aa[i - 1];
  }
  
  cout << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(4);
  cout << endl; 
  /*
  for(int i = 0; i < mp; ++i) {
      for(int j = 0; j < np; ++j) {      
          cout<< am[ i + j * mp] <<' ';
      }
      cout << endl;
  }
  cout << endl;
  */
  int icount = 0;
  do {
     simplx(am, &m, &n, &mp, &np,
            &m1, &m2, &m3, &icase, izrov, iposv,
            &nmax, l1, &mmax, l2, l3);

     if( icase == -1) {
         for(int i = 0; i < mp * np; ++i) am[i] = 0.0;
         am[ 0 + 0 * mp] =  0.0;
         for(int j = 1; j < np; ++j)  am[ 0 + j * mp] =  1.0;
         for(int i = 1; i <= m1; ++i) {
             am[ i + 0 * mp] =  upper_bound;
             am[ i + i * mp] = -1.0;
         }
         am[m1 + 1 + 0 * mp] = res;
         for(int i = 1; i <= n; ++i) 
         am[m1 + 1 + i * mp] = -aa[i - 1];
         /*
         for(int i = 0; i < mp; ++i) {
             for(int j = 0; j < np; ++j) {      
                 cout<< am[ i + j * mp] <<' ';
             }
             cout << endl;
         }
         cout << endl;
         */
	 int scheme = 2;
	 if(scheme == 1) { // increase upper bound
            double upper_bound_new;	     
            double tmp13 = -aa[0] - aa[2];
            double tmp24 = -aa[1] - aa[3];
            if(tmp13 < 0 ) {
	       upper_bound_new = res / fabs(tmp13);
               am[ 1 + 0 * mp] =  upper_bound_new;
	       am[ 3 + 0 * mp] =  upper_bound_new;
	    } else {
	       upper_bound_new = res / fabs(tmp24);
               am[ 2 + 0 * mp] =  upper_bound_new;
	       am[ 4 + 0 * mp] =  upper_bound_new;
	    }
	 } else { // reduce res.
            double tmp13 = -aa[0] - aa[2];
            double tmp24 = -aa[1] - aa[3];
            if(tmp13 < 0 ) {
	       res = fabs(tmp13);
	    } else {
               res = fabs(tmp24);
	    }

	    cout << " scheme2: "
		 << res   << ' '
		 << aa[0] << ' ' << aa[1] << ' '
		 << aa[2] << ' ' << aa[3] << endl;
	 }
	 ++icount;
	 
	 //cout <<"icount = " << icount << ' '
         //     << upper_bound_new << endl;
	 //exit(0);
     }
  } while (icase == -1 && icount < 10);
  
  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
     exit(0);
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
     exit(0);
  } else {
     /*
     for(int i = 0; i < mp; ++i) {
         for(int j = 0; j < np; ++j) {      
             cout<< am[ i + j * mp] <<' ';
         }
         cout << endl;
     }
     cout << endl;
     */
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) {
            x[ iposv[i] - 1 ] = am[ i  + 1];
	 }
     }
     for(int i = 0; i < n; ++i)  {
 	 cout << "x[" << i << "] = "
	      << x[i] << endl;
     }
     cout << endl;
     
     double sum = res;
     for(int i = 0; i < n; ++i) {
	 //cout << sum   << ' '
         //     << aa[i] << ' '
         //     << x [i] << ' '
         //     << aa[i] * x[i]<<endl;
         sum -= aa[i] * x[i];
     } 
     cout <<"check : " <<  sum <<" has to be 0!" << endl;
     //for(int j = 0; j < n; ++j) cout << izrov[j]<<' ';
     cout<< endl;
  }
  //exit(0);
  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

void PVar_Sensitivity::calcPermStdBySimplx(double *amR, double *amO, double *xSolu) {
  int icase;
  int n  = nNodePerm;
  int np = n + 1;   
  int m1 = n;
  int m2 = 0; 
  int m3 = 1;
  int m  = m1 + m2 + m3; 
  int mp = m + 2;
  int nm1m2 = n + m1 + m2;
  int m1p1 = m1 + 1;
  double * am    = new double[mp * np];
  int    * izrov = new int [n];
  int    * iposv = new int [m];

  // Initialize to 0.
  for(int i = 0; i < mp * np; ++i) am[i] = 0.0;
  
  // The first row in the tableau, or row 0
  am[ 0 + 0 * mp] =  0.0;
  for(int j = 1; j < np; ++j) {
      //am[ 0 + j * mp] =  1.0;
      am[ 0 + j * mp] =  amO[j - 1];
  }
 
  // Rows from 1 to m1 in the tableau, the part with <=  
  for(int i = 1; i <= m1; ++i) {
      am[ i + 0 * mp] =  1.0;
      am[ i + i * mp] = -1.0;
  }
  
  for(int i = 0; i < m3; ++i) {
          am[ m1p1 + i + 0 * mp] = 0.0;
      for(int j = 1; j <= np; ++j) {
          am[ m1p1 + i + j * mp] = amR[j - 1];
	  //cout << amR[j - 1] << endl;
      }
  }
  int nmax = n;
  int mmax = m;
  int    * l1    = new int [nmax];
  int    * l2    = new int [mmax];
  int    * l3    = new int [mmax];
  simplx(am, &m, &n, &mp, &np,
         &m1, &m2, &m3, &icase, izrov, iposv,
         &nmax, l1, &mmax, l2, l3);

  if(        icase ==  1) {
     cout << "Unbounded objective function " << endl;
  } else if (icase == -1) {
     cout << "No solution satisfy constraints given " << endl;
  } else {
     for(int i = 0; i < m; ++i) {
         if(iposv[i] <= n) {
            //cout << "x[" << iposv[i] << "] = "<< am[ i  + 1] << endl;
	    xSolu[iposv[i] - 1] = am[ i  + 1];
	 }
     }
     //cout << endl;
  }
  
  /*
  // Using for checking !!
  cout << endl << endl;
  double sum = 0;
  for(int i = 0; i < n; ++i) {
      sum += xSolu[i] * amR[i];
      cout << "i = " << i << ' ' 
	   << xSolu[i]    << ' ' 
	   << amR[i]      << ' '
	   << xSolu[i] * amR[i] << endl;
  }
  cout << "sum = " << sum << endl;
  exit(0);
  */
  
  delete [] izrov;
  delete [] iposv;
  delete [] am;
  delete [] l1;
  delete [] l2;
  delete [] l3;
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{ 
  void simplx_(double *amatrix, int *m, int *n, int *mp, int *np,
               int *m1, int *m2, int *m3, int *icase, int *izrov, int *iposv,
               int *nmax, int *l1, int* mmax, int *l2, int *l3);
};

void PVar_Sensitivity::simplx(double *amatrix, int *m, int *n, int *mp, int *np, int *m1, 
                              int *m2, int *m3, int *icase, int *izrov, int *iposv,
                              int *nmax, int *l1, int* mmax, int *l2, int *l3) 
{
  simplx_(amatrix, m, n, mp, np, m1, m2, m3, icase, izrov, iposv,
          nmax, l1, mmax, l2, l3);
}

