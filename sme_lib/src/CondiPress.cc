#include "CondiPress.h"

/*
CondiPress::CondiPress(Grid* grid, Perm *perm, const P0Eqn &P0ee) 
 : gridPtr(grid), permPtr(perm), P0e( &P0ee )
{
   cout<<" CondiPress::CondiPress()"<<endl;
   p0  = P0e->getP0();
   readData();
   initialization();
}
*/

CondiPress::CondiPress(Grid* grid, Perm *perm, P0Eqn* P0ee, CYPEqn *CYPee, 
                       CPPEqn *CPPee, bool debug_) 
 : gridPtr(grid), permPtr(perm), P0e( P0ee ), CYPe (CYPee), CPPe (CPPee),
   debug( debug_ )
{  
   if(debug) {
      cout<< endl 
          << " CondiPress::CondiPress()"<<endl;
   }
   
   p0  = P0e->getP0();

   readData();
   initialization();
}

CondiPress::~CondiPress(){
   if(debug) 
      cout<<"~CondiPress::CondiPress()"<<endl;

   delete mptsPtr;
   delete krigPtr;
   delete prsmPtr;
   delete p0snPtr;
    
   delLocalArrays();
}

void CondiPress::readData() {
   ifstream in_Pcondi;
   in_Pcondi.open("press_condi.in",ios::in);
   if(in_Pcondi.bad() ) {
      cerr << " Can not open file " <<  "press_condi.in"  << endl;
      exit(8);
   } else {
      if(debug) 
         cout << "File press_condi.in was opened successfully!" <<endl;
   }

   //1) masterpoints object
   long seed;
   int num_i_ms, num_j_ms, num_k_ms;
   Junk(in_Pcondi); in_Pcondi >> seed;
   if(debug) cout<<"seed = "<< seed<<endl;
   
   Junk(in_Pcondi); in_Pcondi >> num_i_ms >> num_j_ms >> num_k_ms;
   
   if(debug) cout<<"num_i_ms = " << num_i_ms << ' '<< num_j_ms << ' '
                                 << num_k_ms << endl;
   num_mpts = num_i_ms * num_j_ms * num_k_ms;

   addMasterPointObj(seed, num_i_ms, num_j_ms, num_k_ms);

   //2) pressure measurements object
   // int num_p_meas;
   int i_pmeas, j_pmeas, k_pmeas;
   double p_meas, p_pred;
   
   Junk(in_Pcondi); in_Pcondi >> num_p_meas;

   if(debug) cout<<"num_p_meas = " << num_p_meas <<endl;
   addPressMeasureObj(num_p_meas);
   for(int m = 0; m < num_p_meas; ++m) {
       Junk(in_Pcondi); 
       in_Pcondi >> i_pmeas>>j_pmeas>>k_pmeas>>p_meas;
       if(debug) cout << i_pmeas<<' '<< j_pmeas<<' '<< k_pmeas<<' '
                 << p_meas << endl;
       prsmPtr->setPdata(m, --i_pmeas, --j_pmeas, --k_pmeas, p_meas);
   }
 
   // maximum number of iterations 
   Junk(in_Pcondi); 
   in_Pcondi >> max_iter >> min_iter >> num_perm_meas >> num_step 
             >> num_iter_control;
   if(debug) cout << max_iter      << ' ' 
                  << min_iter      << ' '
                  << num_perm_meas << ' '
                  << num_step      << endl;
   Junk(in_Pcondi);
   in_Pcondi >>      ind_plot >> ind_fobj >> ind_dbug >> ldbg;
   if(debug) cout << ind_plot << ' ' 
                  << ind_fobj << ' ' 
                  << ind_dbug << ' '
                  << ldbg     << endl;
   Junk(in_Pcondi);
   in_Pcondi >> eps1 >> eps3 >> eps5 >> relax >> dconve >> multplier;
   if(debug) cout << eps1      << ' ' 
                  << eps3      << ' '
                  << eps5      << ' ' 
                  << relax     << ' '
                  << dconve    << ' ' 
                  << multplier << endl;
}

void CondiPress::initialization() {
   addKrigWeighObj();
   addP0SensitivityObj();
   addLocalArrays();
}

// =====   Section 1   =====
// ===== add functions =====
void CondiPress::addMasterPointObj(long seed, int num_i_ms, 
                                   int num_j_ms, int num_k_ms)
{
   mptsPtr = new MasterPoint(seed, num_i_ms, num_j_ms, num_k_ms);
}

void CondiPress::addPressMeasureObj(int num) {
   prsmPtr = new PressMeasure(num, gridPtr, debug);
}


void CondiPress::addKrigWeighObj() {
   krigPtr = new KrigWeigh(gridPtr, permPtr, mptsPtr->getLength());
}

void CondiPress::addP0SensitivityObj() {
   p0snPtr = new P0_Sensitivity( *((P0Eqn*) P0e), mptsPtr, 
                                 prsmPtr, krigPtr );
}

void CondiPress::addLocalArrays() {
   Y_avg_UnCd = new double [ gridPtr->getNumPermNode() ];
   Y_std_UnCd = new double [ gridPtr->getNumPermNode() ];
   
   Y_avg_Mast = new double [ num_mpts ];
   Y_std_Mast = new double [ num_mpts ];
   Y_avg_Cond = new double [ num_mpts ];
   Y_std_Cond = new double [ num_mpts ];
   dY_avg     = new double [ num_mpts ];
   dY_std     = new double [ num_mpts ];

   dpAvg_dY   = new double [ num_mpts * num_p_meas];
   dpStd_dY   = new double [ num_mpts * num_p_meas];
   p_avg_UnCd = new double [ gridPtr->getnNode() ];
   p_std_UnCd = new double [ gridPtr->getnNode() ];
   p_avg_Cond = new double [ gridPtr->getnNode() ];
   p_std_Cond = new double [ gridPtr->getnNode() ];
   p_avg_Meas = new double [ gridPtr->getnNode() ];
   p_std_Meas = new double [ gridPtr->getnNode() ];
   p_avg_Old  = new double [ gridPtr->getnNode() ];
   p_std_Old  = new double [ gridPtr->getnNode() ];
   p_avg_New  = new double [ gridPtr->getnNode() ];
   p_std_New  = new double [ gridPtr->getnNode() ];

   cpq        = new double [ num_mpts * num_mpts ];
   apq        = new double [ num_mpts * 2 ];
   bpq        = new double [ num_mpts * 2 ];
   ppq        = new double [ num_mpts ]; 
  
   cpq_var  = new double [ num_mpts * num_mpts ];
   bpq_var  = new double [ num_mpts * 2 ];
   ppq_var  = new double [ num_mpts ];
}

void CondiPress::delLocalArrays() {
   delete [] Y_avg_UnCd;
   delete [] Y_std_UnCd;
   delete [] Y_avg_Mast;
   delete [] Y_std_Mast;
   delete [] Y_avg_Cond;
   delete [] Y_std_Cond;
   delete [] dY_avg;
   delete [] dY_std;
   delete [] dpAvg_dY;
   delete [] dpStd_dY;
   
   delete [] p_avg_UnCd;
   delete [] p_std_UnCd;
   delete [] p_avg_Cond;
   delete [] p_std_Cond;
   delete [] p_avg_Meas;
   delete [] p_std_Meas;
   delete [] p_avg_Old;
   delete [] p_std_Old;
   delete [] p_avg_New;
   delete [] p_std_New;

   delete [] apq;
   delete [] cpq;
   delete [] bpq;
   delete [] ppq;
   
   delete [] cpq_var;
   delete [] bpq_var;
   delete [] ppq_var;
}

// ===== Section 2 (!!! Conditioning Mean Only !!!) =====
void CondiPress::solve_p0(int nWells, const Well *wellArr, const Control &c, 
                          Solver &s, int time_step ) {
        
   //1) set both Unconditonal Y and P moments in this object!;
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output("UnCondP", p_avg_UnCd, p_std_UnCd);
      
   //2) create MasterPoints and calculate the Krigging weights
   krigAtMasterPoints();

   //4) assing p0 into prsmObj
   set_P0_in_prsmPtr();

   //solve the sensititities of p0 and p_var
   p0snPtr->solve(c, s);
   //  p0snPtr->display();

   sumdiv     = prsmPtr->getAbsDiff(sumwall);
   
   cout << "Averaged pressure deviation: " <<endl;
   cout << "sumdiv     = " << sumdiv     << endl;
   cout << "sumwall    = " << sumwall << endl;                            
   if( sumdiv <= .000001 ) {
       return ;
   }
           
   fobj0     = prsmPtr->getSqdDiff();
   fobj_ini  = fobj0;
   cout << "fobj0     = " << fobj0 << endl;
   cout << "Init obj. function " 
        << fobj0/fobj_ini << ' ' 
        << endl;
   
   if( fobj0 < .001 ) {
      cout << "***** Initial Field Close Enough *****" << endl;
      return;
   }

   for(int iter = 0; iter < max_iter; ++iter) {
       cout << endl;
       cout << "!!! ITERATION = " << iter << " !!!" <<endl;
       cout << endl;
 
       setYCondAtMstPts();

       calcPenaltyF();
       minimization();

       dYWithrelax();
       
       permPtr->PermCondition( num_mpts, krigPtr->getWeigh(),
                               dY_avg, dY_std, mptsPtr->getIJK()
                             );
       // solve P0 solution
       P0e->solve(nWells, wellArr, c, s, 0);

       krigAtMasterPoints();
       set_P0_in_prsmPtr();
       
       //solve the sensititities of p0
       p0snPtr->solve(c, s);
       
       fobj0    = prsmPtr->getSqdDiff();
       cout << "Updated obj. function " 
            << fobj0/fobj_ini << ' '
            << iter << endl;

       if( fobj0/fobj_ini < dconve ) {
           cout << endl;
           cout <<"Obj value is found!!" << endl;
           cout << fobj0/fobj_ini << " < " << dconve << endl;
           return;
       }
   }
}

// ===== setUnCondMoments() =====
void CondiPress::setUnCondMoments() {
   for(int i = 0; i < gridPtr->getNumPermNode(); ++i ) {
       Y_avg_UnCd[i] =       permPtr->getYAvg( i );
       Y_std_UnCd[i] = sqrt( permPtr->getYVar( i ) ); 
   }
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_UnCd[i] =   p0[i];
       p_std_UnCd[i] =  CPPe->getUnCdPStd( i );
   }
}

// ===== setYCondAtMstPts() =====
void CondiPress::setYCondAtMstPts() {
   for(int mst = 0; mst < num_mpts; ++mst ) {
       Y_avg_Cond[mst] =       permPtr->getYAvg( mptsPtr->getIm(mst), 
                                                 mptsPtr->getJm(mst), 
                                                 mptsPtr->getKm(mst) 
                                               );
       Y_std_Cond[mst] = sqrt( permPtr->getYVar( mptsPtr->getIm(mst), 
                                                 mptsPtr->getJm(mst), 
                                                 mptsPtr->getKm(mst) 
                                               )
                              );
       if(debug) 
          cout<< "<Y>c =    " << Y_avg_Cond[mst] << ' '
              << "std_Y_c = " << Y_std_Cond[mst] << endl;
   }
}

// ===== set_P0_in_prsmPtr() =====
void CondiPress::set_P0_in_prsmPtr() {
   for(int m = 0; m < prsmPtr->getLength(); ++m) {
       int ijk = gridPtr->getIndex( prsmPtr->getIpmeas(m), 
                                    prsmPtr->getJpmeas(m),
                                    prsmPtr->getKpmeas(m) 
                                  );
       prsmPtr->setPcond(m, p0[ijk]);
   }
}

// ===== set_P0_in_prsmPtr() =====
void CondiPress::set_P0_in_prsmPtr(double * pavg, double* pstd) {
   for(int m = 0; m < prsmPtr->getLength(); ++m) {
       int ijk = gridPtr->getIndex( prsmPtr->getIpmeas(m), 
                                    prsmPtr->getJpmeas(m),
                                    prsmPtr->getKpmeas(m) 
                                  );
       prsmPtr->setPcond(m, pavg[ijk], pstd[ijk]);
   }
}

// ===== dYWithrelax() =====
void CondiPress::dYWithrelax() {
   for(int mst = 0; mst < num_mpts; ++mst ) {
       dY_avg[mst] *= relax;
       dY_std[mst] *= relax;
   }
}

// ===== Section 3 Solve() =====
void CondiPress::solve(int nWells, const Well *wellArr, const Control &c, 
                       Solver &s, int time_step ) {

   //1) set both Unconditonal Y and P moments in this object!;
   //   and write UnCd Pressure Standard Dev.
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output("UnCondP", p_avg_UnCd, p_std_UnCd);

   //1.1) 
   prsmPtr->output("MeaSurP");
   
   //1.2) calculate CPP conditional to Pressure measurements, we
   //   have two algorithm(), one is linear interpretation,
   //   the other is the fixed pressure
   int algorithm = 2;
   if( algorithm == 1) {
       condCPPbyGaussDist();
   } else {
       condCPPbySolveEqn(nWells, wellArr, c, s);
   }
   output("TarGetP", p_avg_Meas, p_std_Meas);

   //2) create MasterPoints and calculate the Krigging weights
   krigAtMasterPoints();

   //3) assing p0 into prsmObj
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Cond[ i ] = p_avg_UnCd[ i ];
       p_std_Cond[ i ] = p_std_UnCd[ i ];
   }
   set_P0_in_prsmPtr(p_avg_Cond, p_std_Cond);

   //4) solve the sensititities of p0 and p_std
   //P0e->solve(nWells, wellArr, c, s, 0);
   //p0snPtr->solve(c, s);
   calcSensitivity(nWells, wellArr, c, s, p_avg_Cond, p_std_Cond);
   sumdiv     = prsmPtr->getAbsDiff_D(sumwall);

   output("ConditP", p_avg_Cond, p_std_Cond);
   
   cout << "Averaged pressure deviation: " <<endl;
   cout << "sumdiv     = " << sumdiv     << endl;
   cout << "sumwall    = " << sumwall << endl;                            
   if( sumdiv <= .000001 ) {
       return ;
   }
           
   fobj0     = prsmPtr->getSqdDiff();
   fobj_ini  = fobj0;
   cout << "fobj0     = " << fobj0 << endl;
   cout << "Init obj. function " 
        << fobj0/fobj_ini << ' ' 
        << endl;
   
   if( fobj0 < .001 ) {
      cout << "***** Initial Field Close Enough *****" << endl;
      return;
   }
   
   for(int iter = 0; iter < max_iter; ++iter) {
       cout << endl;
       cout << "!!! ITERATION = " << iter << " !!!" <<endl;
       cout << endl;
 
       setYCondAtMstPts();
       calcPenaltyFD();
       minimization();
       dYWithrelax();
       
       permPtr->PermCondition( num_mpts, krigPtr->getWeigh(),
                               dY_avg, dY_std, mptsPtr->getIJK()
                             );

       // solve P0 and PStd solution
        P0e->solve(nWells, wellArr, c, s, 0);
       CYPe->solve(nWells, wellArr, c, s, 0);
       CPPe->solve(nWells, wellArr, c, s, 0);
       CPPe->calcCondPStd();
       for(int i = 0; i < gridPtr->getnNode(); ++i ) {
           p_avg_Cond[ i ] = p0[i];
           p_std_Cond[ i ] = CPPe->getCondPStd( i );
       }
        
       krigAtMasterPoints();
       set_P0_in_prsmPtr(p_avg_Cond, p_std_Cond);
       
       //solve the sensititities of p0 and p_std
       // p0snPtr->solve(c, s);
       calcSensitivity(nWells, wellArr, c, s, p_avg_Cond, p_std_Cond);
 
       fobj0    = prsmPtr->getSqdDiff();

       cout << "Updated obj. function " 
            << fobj0/fobj_ini << ' '
            << iter << endl;

       if( fobj0/fobj_ini < dconve ) {
           cout << endl;
           cout <<"Obj value is found!!" << endl;
           cout << fobj0/fobj_ini << " < " << dconve << endl;
	   output("ConditP", p_avg_Cond, p_std_Cond);
           return;
       }
   }
   output("ConditP", p_avg_Cond, p_std_Cond);
}

// ===== Section 3 Solve() =====
void CondiPress::try_Dagan(int nWells, const Well *wellArr, const Control &c, 
                           Solver &s, int time_step ) {
   cout << "Try Dagan ()" << endl;
   //1) set both Unconditonal Y and P moments in this object!;
   //   and write UnCd Pressure Standard Dev.
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output("UnCondP", p_avg_UnCd, p_std_UnCd);

   //1.1) 
   prsmPtr->output("MeaSurP");
   
   //1.2) krigging
   int msu_ijk = prsmPtr->getIJK(0);
   double pvar = CPPe->getCPP(msu_ijk , msu_ijk);
   double * CYY = new double[ gridPtr->getNumPermNode() * gridPtr->getNumPermNode()];
   double * rho = new double[ gridPtr->getNumPermNode() * gridPtr->getnNode()]; 
   double * Yavg= new double[ gridPtr->getNumPermNode() ];
   for(int k2 = 0; k2 < gridPtr->getPermNz(); ++k2) {
       for(int j2 = 0; j2 < gridPtr->getPermNy(); ++j2) {
           for(int i2 = 0; i2 < gridPtr->getPermNx(); ++i2) {
	       int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
               for(int k1 = 0; k1 < gridPtr->getPermNz(); ++k1) {
                   for(int j1 = 0; j1 < gridPtr->getPermNy(); ++j1) {
                       for(int i1 = 0; i1 < gridPtr->getPermNx(); ++i1) {
			   int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
                           CYY[ijk1 + ijk2 * gridPtr->getNumPermNode()]
			   = permPtr->getCYY(i1, j1, k1, i2, j2, k2) ;
		       }     
		   }
	       }
               rho[ ijk2 ] = CYPe->getCYP(msu_ijk, ijk2)/ pvar;
	   }
       }
   }

   double diff = prsmPtr->getPmeas(0) - p0[msu_ijk];
   for(int k2 = 0; k2 < gridPtr->getPermNz(); ++k2) {
       for(int j2 = 0; j2 < gridPtr->getPermNy(); ++j2) {
           for(int i2 = 0; i2 < gridPtr->getPermNx(); ++i2) {
	       int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
	       Yavg[ijk2] = permPtr->getYAvg(ijk2) + rho[ ijk2 ] * diff;
	       permPtr->setYAvg(ijk2, Yavg[ijk2]);
 	   }
       }
   }
   
   for(int k2 = 0; k2 < gridPtr->getPermNz(); ++k2) {
       for(int j2 = 0; j2 < gridPtr->getPermNy(); ++j2) {
           for(int i2 = 0; i2 < gridPtr->getPermNx(); ++i2) {
	       int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
               for(int k1 = 0; k1 < gridPtr->getPermNz(); ++k1) {
                   for(int j1 = 0; j1 < gridPtr->getPermNy(); ++j1) {
                       for(int i1 = 0; i1 < gridPtr->getPermNx(); ++i1) {
			   int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
			   double value = rho[ ijk2 ] * CYPe->getCYP(msu_ijk, ijk1);
			   permPtr->setCYY_Cond(i1, j1, k1, i2, j2, k2, value);
                           //CYY[ijk1 + ijk2 * gridPtr->getNumPermNode()]
			   // -= rho[ ijk2 ] * CYPe->getCYP(msu_ijk, ijk1);  
		       }     
		   }
	       }
	   }
       }
   }
    
    P0e->solve(nWells, wellArr, c, s, 0);
   CYPe->solve(nWells, wellArr, c, s, 0);
   CPPe->solve(nWells, wellArr, c, s, 0);
   CPPe->calcCondPStd();
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Cond[ i ] = p0[i];
       p_std_Cond[ i ] = CPPe->getCondPStd( i );
   }
   output("ConditP", p_avg_Cond, p_std_Cond);

   delete [] rho;
   delete [] CYY;
   delete [] Yavg;
   
   //1.3) calculate CPP conditional to Pressure measurements, we
   //   have two algorithm(), one is linear interpretation,
   //   the other is the fixed pressure
   int algorithm = 2;
   if( algorithm == 1) {
       condCPPbyGaussDist();
   } else {
       condCPPbySolveEqn(nWells, wellArr, c, s);
   }
   output("TarGetP", p_avg_Meas, p_std_Meas);
   
}

/*
// ===== Section 3 Solve() =====
void CondiPress::solve(int nWells, const Well *wellArr, const Control &c, 
                       Solver &s, int time_step ) {

   //1) set both Unconditonal Y and P moments in this object!;
   //   and write UnCd Pressure Standard Dev.
   CPPe->calcUnCdPStd();
   setUnCondMoments();
   output("UnCondP", p_avg_UnCd, p_std_UnCd);

   //1.1) 
   prsmPtr->output("MeaSurP");
   
   //1.2) calculate CPP conditional to Pressure measurements, we
   //   have two algorithm(), one is linear interpretation,
   //   the other is the fixed pressure
   int algorithm = 2;
   if( algorithm == 1) {
       condCPPbyGaussDist();
   } else {
       condCPPbySolveEqn(nWells, wellArr, c, s);
   }
   output("TarGetP", p_avg_Meas, p_std_Meas);

   //2) create MasterPoints and calculate the Krigging weights
   krigAtMasterPoints();

   //3) assing p0 into prsmObj
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Cond[ i ] = p_avg_UnCd[i];
       p_std_Cond[ i ] = p_std_UnCd[i];
       //p_std_Cond[ i ] = 0.;
   }
   set_P0_in_prsmPtr(p_avg_Cond, p_std_Cond);
   //set_P0_in_prsmPtr(p_avg_UnCd, p_std_UnCd);

   //4) solve the sensititities of p0 and p_std
   P0e->solve(nWells, wellArr, c, s, 0);
   p0snPtr->solve(c, s);
   calcSensitivity(nWells, wellArr, c, s, p_avg_Cond, p_std_Cond);
   sumdiv     = prsmPtr->getAbsDiff_D(sumwall);

   output("ConditP", p_avg_Cond, p_std_Cond);
   
   cout << "Averaged pressure deviation: " <<endl;
   cout << "sumdiv     = " << sumdiv     << endl;
   cout << "sumwall    = " << sumwall << endl;                            
   if( sumdiv <= .000001 ) {
       return ;
   }
           
   fobj0     = prsmPtr->getSqdDiff();
   fobj_ini  = fobj0;
   cout << "fobj0     = " << fobj0 << endl;
   cout << "Init obj. function " 
        << fobj0/fobj_ini << ' ' 
        << endl;
   
   if( fobj0 < .001 ) {
      cout << "***** Initial Field Close Enough *****" << endl;
      return;
   }
   
   for(int iter = 0; iter < max_iter; ++iter) {
       cout << endl;
       cout << "!!! ITERATION = " << iter << " !!!" <<endl;
       cout << endl;
 
       setYCondAtMstPts();
       calcPenaltyFD();
       minimization();
        
       dYWithrelax();
       
       permPtr->PermCondition( num_mpts, krigPtr->getWeigh(),
                               dY_avg, dY_std, mptsPtr->getIJK()
                             );

       // solve P0 and PStd solution
        P0e->solve(nWells, wellArr, c, s, 0);
       CYPe->solve(nWells, wellArr, c, s, 0);
       CPPe->solve(nWells, wellArr, c, s, 0);
       CPPe->calcCondPStd();
       for(int i = 0; i < gridPtr->getnNode(); ++i ) {
           p_avg_Cond[ i ] = p0[i];
           p_std_Cond[ i ] = CPPe->getCondPStd( i );
       }
        
       krigAtMasterPoints();
       set_P0_in_prsmPtr(p_avg_Cond, p_std_Cond);
       
       //solve the sensititities of p0 and p_std
       p0snPtr->solve(c, s);
       calcSensitivity(nWells, wellArr, c, s, p_avg_Cond, p_std_Cond);
 
       fobj0    = prsmPtr->getSqdDiff();

       cout << "Updated obj. function " 
            << fobj0/fobj_ini << ' '
            << iter << endl;

       if( fobj0/fobj_ini < dconve ) {
           cout << endl;
           cout <<"Obj value is found!!" << endl;
           cout << fobj0/fobj_ini << " < " << dconve << endl;
	   output("ConditP", p_avg_Cond, p_std_Cond);
           return;
       }
   }
   output("ConditP", p_avg_Cond, p_std_Cond);
}
*/

void CondiPress::calcSensitivity(int nWells, const Well *w, const Control &c,
                                 Solver &s, double* pavg, double *pstd) {
                                 
   if(debug) cout << "CalcSensitivity( )" << endl;
   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Old[ i ] = pavg[i];
       p_std_Old[ i ] = pstd[i];
   }

   //output("Pre", p_avg_Old, p_std_Old);
   double dY_avg_tmp = 0.05;
   for(int mst = 0; mst < mptsPtr->getLength(); mst++) {
       //cout << mptsPtr->getIm(mst)<<' ' << mptsPtr->getJm(mst)<<' '
       //     << mptsPtr->getKm(mst)<<' ' << mptsPtr->getIJK(mst)<< endl;

       //1) copy to tmp
       permPtr->Copy_PermYAvgToTmp();
       
       //2) set Perm
       permPtr->SetPermYAvg( mst, num_mpts, 
                             dY_avg_tmp,krigPtr->getWeigh() 
                           );
       
       //permPtr->wrtPermMomnt();
       
       //3) simulation
        P0e->solve(nWells, w, c, s, 0);
       CYPe->solve(nWells, w, c, s, 0);
       CPPe->solve(nWells, w, c, s, 0);

       //CPPe->calcCondPStd();
       CPPe->calcUnCdPStd();        
       for(int i = 0; i < gridPtr->getnNode(); ++i ) {
           p_avg_New[ i ] = p0[i];
           p_std_New[ i ] = CPPe->getUnCdPStd( i );
       }
       
       //output("Post", p_avg_New, p_std_New);
       for(int i = 0; i < num_p_meas; ++i) {
           int ijk = gridPtr->getIndex( prsmPtr->getIpmeas(i),
                                           prsmPtr->getJpmeas(i),
                                        prsmPtr->getKpmeas(i)
                                      );
           dpAvg_dY[ i + mst * num_p_meas] = (p_avg_New[ ijk ] - p_avg_Old[ ijk ])
                                           / dY_avg_tmp;
           dpStd_dY[ i + mst * num_p_meas] = (p_std_New[ ijk ] - p_std_Old[ ijk ])
                                           / dY_avg_tmp;
           cout << "m = " << mst << ' ' << dpAvg_dY[ i + mst * num_p_meas] << ' '
                << p0snPtr->getP0Sensitivity(i, mst) << ' '
                << dpStd_dY[ i + mst * num_p_meas] << endl;
       }
         
       //4) copy tmp back
       permPtr->Copy_TmpToPermYAvg();
   }
   
}

void CondiPress::set_PStd_in_pstdPtr() {
   if(debug) cout<<"AFTER Condition or set_PStd_in_pstdPtr()"
                 <<endl;
   for(int m = 0; m < pstdPtr->getLength(); ++m) {
       int ijk = gridPtr->getIndex( pstdPtr->getIpmeas(m), 
                                    pstdPtr->getJpmeas(m),
                                    pstdPtr->getKpmeas(m) 
                                  );
       
       if(debug)
          cout << pstdPtr->getIpmeas(m) << ' '
               << pstdPtr->getJpmeas(m) << ' '
               << pstdPtr->getKpmeas(m) << ' '
               << pVsnPtr->getPStd(ijk) << endl;
   
       pstdPtr->setPcond(m, pVsnPtr->getPStd(ijk));
   }
}

void CondiPress::condCPPbySolveEqn(int nWells, const Well *w, 
                                   const Control &c, Solver &s) {
   // condition part
   P0e->solve(nWells, w, c, s, 0, num_p_meas, 
              prsmPtr->getIpmeas(),
              prsmPtr->getJpmeas(),
              prsmPtr->getKpmeas(),
              prsmPtr->getPmeas()
             );
   CYPe->solve(nWells, w, c, s, 0);
   CPPe->solve(nWells, w, c, s, 0);
   CPPe->calcCondPStd();

   for(int i = 0; i < gridPtr->getnNode(); ++i ) {
       p_avg_Meas[i] = p0[i];
       p_std_Meas[i] = CPPe->getCondPStd( i );
   }
}

void CondiPress::condCPPbyGaussDist() {
   krigPtr->condGaussDist(gridPtr->getnNode(), num_p_meas, 
                          prsmPtr->getIJKpmeas(), CPPe->getCPP(),
                          p_std_Meas
                         );
}

void CondiPress::krigAtMasterPoints() {
   
   mptsPtr->calMasterPts( gridPtr->getPermNx(), 
                          gridPtr->getPermNy(),
                          gridPtr->getPermNz(),
                          gridPtr->getPermX(), 
                          gridPtr->getPermY(),
                          gridPtr->getPermZ() 
                        );
   setYUnCdAtMstPts();
   
   krigPtr->calWeigh( mptsPtr->getIm(),
                      mptsPtr->getJm(), 
                      mptsPtr->getKm() 
                    );
   
   /* //try area
   mptsPtr->setMstPts( 
                      gridPtr->getPermX(), 
                      gridPtr->getPermY(),
                      gridPtr->getPermZ()
                     );
   setYUnCdAtMstPts();
   
   krigPtr->calWeigh( mptsPtr->getIm(),
                      mptsPtr->getJm(), 
                      mptsPtr->getKm() 
                    );
   */
}

void CondiPress::setYUnCdAtMstPts() {
   for(int mst = 0; mst < num_mpts; ++mst ) {
       int ijk = gridPtr->getPermIndex(mptsPtr->getIm(mst), 
                                       mptsPtr->getJm(mst), 
                                       mptsPtr->getKm(mst)
                                      );
       Y_avg_Mast[mst] = Y_avg_UnCd[ijk];
       Y_std_Mast[mst] = Y_std_UnCd[ijk];
   }
}

void CondiPress::addPStdSensitivityObj(int num_pStd) {

   pstdPtr = new PressMeasure(num_pStd, gridPtr, debug);
   
   for(int m = 0; m < num_pStd; ++m) {
       int i   = m;
       int j   = m;
       int k   = 0;
       int ijk = gridPtr->getIndex(i, j, k);
       pstdPtr->setPdata(m, i, j, k, p_std_Meas[ijk]);
       //cout << p_std_Meas[ijk] << endl;
   }
   //exit(0);  
   pVsnPtr = new PVar_Sensitivity( *((P0Eqn*) P0e), CYPe, CPPe,
                                   mptsPtr, pstdPtr, krigPtr );
   set_PStd_in_pstdPtr();
   if(debug) pstdPtr->display();
}

void CondiPress::delPStdSensitivityObj() {
   delete pstdPtr;
   delete pVsnPtr;   
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{
  void buildf_(int *mxwell, int *mxtpe, int *mxstep, int *ntmed, int *nmp,
               int *nwell,  int *nstep, double *htj, double *wpms, double *sumdiv,
               double *pms, double *pcal, double *am_y, double *Y_avg_Mast, 
               double *Y_std_Mast, double *cpq, double *ppq, double *apq, 
               double *bpq,double *Y_avg_Cond, double *dY_avg);
  
  void buildfd_(int *mxwell, int *mxtpe, int *mxstep, int *ntmed, int *nmp,
               int *nwell,  int *nstep, double *htj, double *htj1, double *wpms, double *sumdiv,
               double *pms, double *pcal, double *pstd, double *am_y, double *Y_avg_Mast, 
               double *Y_std_Mast, double *cpq, double *ppq, double *apq, 
               double *bpq,double *Y_avg_Cond, double *dY_avg);

  void buildf_var_(int *mxwell, int *mxtpe, int *mxstep, int *ntmed, int *nmp,
               int *nwell,  int *nstep, double *htj, double *wpms, double *sumdiv_var,
               double *pms, double *pcal, double *am_y, double *Y_avg_Mast, 
               double *Y_std_Mast, double *cpq_var, double *ppq_var, double *apq, 
               double *bpq_var, double *Y_std_Cond, double *dY_std);
  
  void minim2_(int *nkrig, int *it_min, int *ifobj, int *idbg, int *ldbg,
               int *mxtpe, double *fobj0, double *sumwall, double *eps1,
               double *eps3, double *eps5, double *cpq, double *ppq,
               double *apq, double *bpq, double *dY_avg);
};

void CondiPress::calcPenaltyFD() {
  buildfd_(&num_p_meas, &num_mpts, &num_step, &num_perm_meas, &num_mpts,
           &num_p_meas, &num_step, dpAvg_dY, dpStd_dY,
           prsmPtr->getWpmeas(), &sumdiv,
           prsmPtr->getPmeas(), prsmPtr->getPpred(), prsmPtr->getPpstd(),
	   &multplier,
           Y_avg_Mast, Y_std_Mast, cpq, ppq, apq, bpq, Y_avg_Cond, dY_avg);
}

void CondiPress::calcPenaltyF() {
  buildf_(&num_p_meas, &num_mpts, &num_step, &num_perm_meas, &num_mpts,
          &num_p_meas, &num_step, p0snPtr->getP0Sensitivity(), 
          prsmPtr->getWpmeas(), &sumdiv,
          prsmPtr->getPmeas(), prsmPtr->getPpred(), &multplier, 
          Y_avg_Mast, Y_std_Mast, cpq, ppq, apq, bpq, Y_avg_Cond, dY_avg);
}

void CondiPress::calcPenaltyFVar() {
  //cout << "calcPenaltyFVar()" << endl;
  //for(int i = 0; i < num_mpts; ++i) {
  //    cout << Y_std_Mast[i] << ' '<< Y_std_Cond[i] << endl;
  // }
  buildf_var_(&num_pstd_meas, &num_mpts, &num_step, &num_perm_meas, &num_mpts,
              &num_pstd_meas, &num_step, pVsnPtr->getPVarSensitivity(),
               pstdPtr->getWpmeas(), &sumdiv_var,
               pstdPtr->getPmeas(),pstdPtr->getPpred(), &multplier,
               Y_avg_Mast, Y_std_Mast,
               cpq_var, ppq_var, apq, bpq_var, Y_std_Cond, dY_std);
}

void CondiPress::minimization() {
  /*
  cout << num_mpts << ' ' << min_iter<< ' '
       << ind_fobj << ' ' << ind_dbug<< ' '
       << ldbg     << ' ' << fobj0   << ' ' 
       << sumwall  << ' ' << eps1    << ' '
       << eps3     << ' ' << eps5    << endl;
  */
  int dummy = ind_fobj;
  minim2_(&num_mpts, &min_iter, &dummy, &ind_dbug, &ldbg, 
          &num_mpts, &fobj0, &sumwall,
          &eps1, &eps3, &eps5, cpq, ppq, apq, bpq, dY_avg);
}

void CondiPress::minimizationVar() {
   int dummy = ind_fobj;
   /*
   for(int mst = 0; mst < num_mpts; ++mst ) 
       cout << dY_std[mst]<<' ';
   cout<<endl;
   */
   minim2_(&num_mpts, &min_iter, &dummy, &ind_dbug, &ldbg, 
           &num_mpts, &fobj0_var, &sumwall,&eps1, &eps3, &eps5,
           cpq_var, ppq_var, apq, bpq_var, dY_std);
   /*
   for(int mst = 0; mst < num_mpts; ++mst ) 
       cout << dY_std[mst]<<' ';
   cout<<endl;
   */
}

//===== output() =====
void CondiPress::output(char *file, double *avg, double *std){
   char * ending1 = "_1d.out";
   char * ending2 = "_2d.out";
   int length  = strlen(file);
   int length1 = strlen(ending1);
   int length2 = strlen(ending2);
   char *file_1d = new char [length + length1];   
   char *file_2d = new char [length + length2];
   strcpy(file_1d, file);
   strcpy(file_2d, file);
   strcat(file_1d,ending1);
   strcat(file_2d,ending2); 
   
   ofstream os1(file_1d, ios::out);
   ofstream os2(file_2d, ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);
   os2 << setiosflags(ios::fixed | ios::showpoint) << setprecision(5);    

   for(int k = 0; k < gridPtr->getNz(); k++) {
       for(int j = 0; j < gridPtr->getNy(); j++) {
           int i = j;
           int ijk = gridPtr->getIndex(i, j, k);
           double dist = sqrt(  gridPtr->getX(i) * gridPtr->getX(i) 
                              + gridPtr->getY(i) * gridPtr->getY(i) );
           os1 << dist << ' '<< avg[ijk] << ' ' << std[ijk] << endl;
               
           for(int i = 0; i < gridPtr->getNx(); i++) {
               int ijk = gridPtr->getIndex(i, j, k);
               os2 << gridPtr->getX(i) << ' ' << gridPtr->getY(j) << ' '
                   << avg[ijk] << ' ' << std[ijk] << endl;
           }
           os2<<endl;
       }
   }
   os1<<endl;
   os2<<endl;
   
   os1.close();
   os2.close();

   delete [] file_1d;
   delete [] file_2d;
}

void CondiPress::output_p0() {
   mptsPtr->output();
   prsmPtr->output();
   krigPtr->output();
   p0snPtr->output();
   
   output("CondiP0", (double*) p0, (double*) p0);
}

void CondiPress::output(){
   mptsPtr->output();
   prsmPtr->output();
   krigPtr->output();
   p0snPtr->output();
   //pnltPtr->output();
   
   ofstream os("CondiPressure_2d.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);
   for(int k = 0; k < gridPtr->getNz(); k++) {
       for(int j = 0; j < gridPtr->getNy(); j++) {
           for(int i = 0; i < gridPtr->getNx(); i++) {
               int ijk = gridPtr->getIndex(i, j, k);
               os << gridPtr->getX(i) << ' ' 
                  << gridPtr->getY(j) << ' '
                  << p0[ijk] << ' '
                  << pVsnPtr->getPStd(ijk)
                  << endl;
           }
           os<<endl;
       } 
   }          
   os<<endl;
   os.close();
   
   ofstream os1("CondiPressure_1d.out", ios::out);
   os1 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(5);
   cout << "print CondiPressure_1d.out is invoked!" << endl;
   for(int k = 0; k < gridPtr->getNz(); k++) {
       for(int j = 0; j < gridPtr->getNy(); j++) {
           int i = j;
           int ijk = gridPtr->getIndex(i, j, k);
           double dist = sqrt(  gridPtr->getX(i) * gridPtr->getX(i) 
                              + gridPtr->getY(i) * gridPtr->getY(i) );
           os1 << dist   << ' '
               << p0[ijk]<< ' '
               << pVsnPtr->getPStd(ijk) 
               << endl;
       } 
   }          
   os1<<endl;
   os1.close();
   
   ofstream os2("CondiPerm_2d.out", ios::out);
   os2 << setiosflags(ios::fixed | ios::showpoint)
       << setprecision(5);
   for(int k = 0; k < gridPtr->getPermNz(); k++) {
       for(int j = 0; j < gridPtr->getPermNy(); j++) {
           for(int i = 0; i < gridPtr->getPermNx(); i++) {
               int ijk = gridPtr->getPermIndex(i, j, k);
               os2 << i << ' ' 
                   << j << ' '
                   << permPtr->getYAvg(i, j, k) <<' '
                   << permPtr->getYVar(i, j, k) <<' '
                   << endl;
           }
           os2<<endl;
       } 
   }          
   os2<<endl;
   os2.close();
}

void CondiPress::outputCondiPerm(){
   ofstream os2("CondiPerm_2d.out", ios::out);
   os2 << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(5);
   for(int k = 0; k < gridPtr->getPermNz(); k++) {
       for(int j = 0; j < gridPtr->getPermNy(); j++) {
           for(int i = 0; i < gridPtr->getPermNx(); i++) {
               int ijk = gridPtr->getPermIndex(i, j, k);
               os2 << i << ' '
                   << j << ' '
                   << permPtr->getYAvg(i, j, k) <<' '
                   << permPtr->getYVar(i, j, k) <<' '
                   << endl;
           }
           os2<<endl;
       }
   }
   os2<<endl;
   os2.close();
}

/*
void CondiPress::addPenaltyObj() {
   pnltPtr = new Penalty(permPtr, mptsPtr, prsmPtr, p0snPtr);
}
*/

