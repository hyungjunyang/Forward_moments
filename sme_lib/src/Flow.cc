#include "Flow.h"

Flow::Flow(char *Dir, Domain *domain, bool debug_) 
:domainPtr(domain), debug(debug_)
{
  if(debug) cout<<"Flow::Flow()" << endl;

  unit    = domainPtr->getUnit();
  gridPtr = domainPtr->getGrid();
  regnPtr = domainPtr->getRegn();
  fluidPtr= domainPtr->getFluid();
  
  nDt     = domainPtr->getDtNum(); 
  contPtr = domainPtr->getContr();

  nWells  = domainPtr->getWellNum();
  wellArr = domainPtr->getWellArr();
  
  initialize( Dir);
}

Flow::~Flow() {
  if(debug) cout<<"Flow::~Flow()" << endl;
  
  delete permPtr;
  delete P2e;
  delete CPPe; 
  delete CYPe; 
  delete P0e;
  delete veloPtr;
  if(ss!=NULL){delete ss;}

}

void Flow::initialize(char *Dir) {
  double startT, endT, T0;
  T0 = startT = clock();
  
  //readUnit  ( is );
  //addGridObj( is , Dir );
  //addRegnObj( is , unit );
  //addFluidObj(is , unit );
  //addWellArr( is );
  //addContObj( is );
  
  addPermObj( unit, Dir );
  //permPtr->readCondYMomntFile(); //Comment out if no conditional perm with perm measurement.
  addPresObj(    );
  addSolverObj(  );
    
  endT = clock();
  cout << " " << endl;
  cerr << "===== Initial Time: " 
       << (endT-startT)/CLOCKS_PER_SEC 
       << " sec =====" << endl;
}

void Flow::addPermObj( Unit unit, char *Dir) {
	permPtr = new Perm(unit, regnPtr, gridPtr, Dir, (contPtr->getFlagPressCond() == 1)||(contPtr->getFlagSatCond() == 1)||(contPtr->getFlagProdCond() == 1), debug);
}


void Flow::addPresObj() {
  P0e  = new P0Eqn(*gridPtr, regnPtr, permPtr, fluidPtr, nWells, 
                   wellArr, *contPtr, unit, debug);
  CYPe = new CYPEqn(*((P0Eqn*) P0e), debug);
  CPPe = new CPPEqn(*((P0Eqn*) P0e), *((CYPEqn*) CYPe), debug);
  P2e  = new P2Eqn (*((P0Eqn*) P0e), *((CYPEqn*) CYPe), debug);
}

void Flow::addSolverObj() {
  cout << endl;
  cout << "===== Solver Info =====" << endl;
  switch( contPtr->getSolverChoice() ) {
  case FULLLU:
    ss = new FullLUSolver( *P0e );
    cout << " Full LU Solver has been created" << endl;
    break;
  case LAPACK_BAND:
    ss = new LapackBandSolver( *P0e );
    cout << " LAPACK Banded Solver has been created" << endl;
    break;
  case LAPACK_FULL:
    ss = new LapackFullSolver( *P0e );
    cout << " LAPACK Full Solver has been created" << endl;
    break;
  case LAPACK_SYM_BAND:
    ss = new LapackSymBandSolver( *P0e );
    cout << " LAPACK Symmetric Banded Solver has been created" << endl;
    break;
  default:
    cerr << "WRONG Solver choice -----\n";
  }
}

void Flow::delSolverObj() {
  delete ss;
  ss = NULL;
}

void Flow::solve(char * Dir, ifstream &is) {

  double startT, endT, T0;
  T0 = startT = clock();

  // --- start calculation ------------------- 
  for(int iter = 0; iter < nDt; iter++) {
      cout << " " << endl;
      cout << "   ===== Time Step " << iter + 1 << " ====="  << endl;
               
      // --- calculate P0 ---
      startT = clock();
      P0e->solve(nWells, wellArr, *contPtr, *ss, iter);
      endT = clock();
      cerr << "     H0  time: " << (endT-startT)/CLOCKS_PER_SEC << " sec"
           << endl;

      //* 
      // --- calculate CYP ---
      startT = clock();
      CYPe->solve(nWells, wellArr, *contPtr, *ss, iter);
      endT = clock();
      cerr << "     CYP time: " << (endT-startT)/CLOCKS_PER_SEC << " sec"
           << endl;
               
      // --- calculate CPP ---
      startT = clock();
      CPPe->solve(nWells, wellArr, *contPtr, *ss, iter);
      //((CPPEqn*) CPPe)->printCPP();
      endT = clock();
      cerr << "     CPP time: " << (endT-startT)/CLOCKS_PER_SEC << " sec"
           << endl;

      // --- calculate P2 ---
      //startT = clock();
      //P2e->solve(nWells, wellArr, *contPtr, *ss, iter);
      //endT = clock();
      //cerr << "     H2time: " << (endT-startT)/CLOCKS_PER_SEC << " sec"
      //    << endl;
  }          
       
  cerr << "     Total Running Time (without Q): " 
       << (endT-T0)/CLOCKS_PER_SEC << " sec" << endl;

  writeToFile(Dir);
  
  // --- now calculate the velocity moments -----
  addVeloObj( Dir );
}

void Flow::addVeloObj(char * Dir) {
  double startT, endT;
  veloPtr = new Velocity(gridPtr, regnPtr, permPtr, fluidPtr,  *((P0Eqn*) P0e), 
                       *((CYPEqn*) CYPe), *((CPPEqn*) CPPe), *((P2Eqn*) P2e), 
                         unit, Dir, debug
                      );

  startT = clock();
  veloPtr->solve( *contPtr, nDt - 1);  
  veloPtr->writeCv(*contPtr);

  endT = clock();
  cerr << endl;
  cerr << "===== Velocity time: " 
       << (endT-startT)/CLOCKS_PER_SEC 
       << " sec =====" << endl;
}

void Flow::writeToFile(char *Dir) {
  // --- output P0, CYP, CPP, P2 ---
  char file[100];
  int slice_index_x, slice_index_y, slice_index_z;
  int flag;

  // (x,y) two dimensional surface plot.
  if( gridPtr->getNx() > 1 && gridPtr->getNy() > 1 && gridPtr->getNz() == 1 ) {
      flag = 0;
      sprintf( file, "%sP0_2D_XY.out", Dir );
      P0e->output( file, flag, *contPtr );
              
      sprintf( file, "%sP2_2D_XY.out", Dir );
      P2e->output( file, flag, *contPtr );
                 
      sprintf( file, "%sCYP_2D_XY.out", Dir );
      CYPe->output( file, flag, *contPtr );

      sprintf( file, "%sCPY_2D_XY.out", Dir );
      ((CYPEqn*)CYPe)->outputCPY( file, flag, *contPtr );
                  
      sprintf( file, "%sCPP_2D_XY.out", Dir );
      CPPe->output( file, flag, *contPtr );
             
      sprintf( file, "%sP_Var_2D_XY.out", Dir );
      ((CPPEqn*) CPPe)->output_PVari( file, flag );
  }        

  // (x,z) two dimensional surface plot.
  if( gridPtr->getNx() > 1 && gridPtr->getNz() > 1 && gridPtr->getNy() == 1 ) {
      flag = 1;
      sprintf( file, "%sP0_2D_XZ.out", Dir );
      P0e->output( file, flag, *contPtr );
      
      sprintf( file, "%sP2_2D_XZ.out", Dir );
      P2e->output( file, flag, *contPtr );
      
      sprintf( file, "%sCYP_2D_XZ.out", Dir );
      CYPe->output( file, flag, *contPtr );
      
      sprintf( file, "%sCPP_2D_XZ.out", Dir );
      CPPe->output( file, flag, *contPtr );

      sprintf( file, "%sP_Var_2D_XZ.out", Dir );
      ((CPPEqn*) CPPe)->output_PVari( file, flag );
  }

  // Diagonal surface plot.
  if( gridPtr->getNx() == gridPtr->getNy() && gridPtr->getNz() > 1 ) {
      flag = 2;
      sprintf( file, "%sP0_2D_XY_Z.out", Dir );
      P0e->output( file, flag, *contPtr );
      
      sprintf( file, "%sP2_2D_XY_Z.out", Dir );
      P2e->output( file, flag, *contPtr );
      
      sprintf( file, "%sCYP_2D_XY_Z.out", Dir );
      CYPe->output( file, flag, *contPtr );
      
      sprintf( file, "%sCPP_2D_XY_Z.out", Dir );
      CPPe->output( file, flag, *contPtr );

      sprintf( file, "%sP_Var_2D_XY_Z.out", Dir );
      ((CPPEqn*) CPPe)->output_PVari( file, flag );
  }

  // Diagonal one dimensional plot coming from (x,y) problem.
  if( gridPtr->getNx() == gridPtr->getNy() && gridPtr->getNz() == 1 ) {
      flag = 3;
      sprintf( file, "%sP0_1D_XY.out", Dir );
      P0e->output( file, flag, *contPtr );
      
      sprintf( file, "%sP2_1D_XY.out", Dir );
      P2e->output( file, flag, *contPtr );
      
      sprintf( file, "%sCYP_1D_XY.out", Dir );
      CYPe->output( file, flag, *contPtr );
      
      sprintf( file, "%sCPP_1D_XY.out", Dir );
      CPPe->output( file, flag, *contPtr );

      sprintf( file, "%sP_Var_1D_XY.out", Dir );
      ((CPPEqn*) CPPe)->output_PVari( file, flag );
  }

  // Diagonal one dimensional plot coming from (x,z) problem.
  if( gridPtr->getNx() == gridPtr->getNz() && gridPtr->getNy() == 1 ) {
      flag = 4;
      sprintf( file, "%sP0_1D_XZ.out", Dir );
      P0e->output( file, flag, *contPtr );
      
      sprintf( file, "%sP2_1D_XZ.out", Dir );
      P2e->output( file, flag, *contPtr );
      
      sprintf( file, "%sCYP_1D_XZ.out", Dir );
      CYPe->output( file, flag, *contPtr );
      
      sprintf( file, "%sCPP_1D_XZ.out", Dir );
      CPPe->output( file, flag, *contPtr );

      sprintf( file, "%sP_Var_1D_XZ.out", Dir );
      ((CPPEqn*) CPPe)->output_PVari( file, flag );
  }

  // Surface plot according to predetermined slice index.
  if( gridPtr->getNx() > 1 && gridPtr->getNy() > 1 && gridPtr->getNz() > 1 ) {
      flag = 9;
      sprintf( file, "%sP0_3D_XYZ.out", Dir );
      P0e->output( file, flag, *contPtr);

      /*	  
      cout << "Need slice index for surface plot ( < 0 for no plot )"
           << endl;
      cout << "Enter slice index in x-dir:... ";
      cin >> slice_index_x;
      if( slice_index_x >= 0 && slice_index_x < gridPtr->getNx() ) {
          flag = 5;
          sprintf( file, "%sP0_2D_YZ_i=%d.out", Dir, slice_index_x );
          P0e->output( file, flag, *contPtr, slice_index_x );
          
          sprintf( file, "%sP2_2D_YZ_i=%d.out", Dir, slice_index_x );
          P2e->output( file, flag, *contPtr, slice_index_x );
      
          sprintf( file, "%sCYP_2D_YZ_i=%d.out", Dir, slice_index_x );
          CYPe->output( file, flag, *contPtr, slice_index_x );
      
          sprintf( file, "%sCPP_2D_YZ_i=%d.out", Dir, slice_index_x );
          CPPe->output( file, flag, *contPtr, slice_index_x );

          sprintf( file, "%sP_Var_2D_YZ_i=%d.out", Dir );
          ((CPPEqn*) CPPe)->output_PVari( file, flag, slice_index_x );
      }

      cout << "Enter slice index in y-dir:... ";
      cin >> slice_index_y;
      if( slice_index_y >= 0 && slice_index_y < gridPtr->getNy() ) {
          flag = 6;

          sprintf( file, "%sP0_2D_XZ_j=%d.out", Dir, slice_index_y );
          P0e->output( file, flag, *contPtr, slice_index_y );
          
          sprintf( file, "%sP2_2D_XZ_j=%d.out", Dir, slice_index_y );
          P2e->output( file, flag, *contPtr, slice_index_y );
      
          sprintf( file, "%sCYP_2D_XZ_j=%d.out", Dir, slice_index_y );
          CYPe->output( file, flag, *contPtr, slice_index_y );
      
          sprintf( file, "%sCPP_2D_XZ_j=%d.out", Dir, slice_index_y );
          CPPe->output( file, flag, *contPtr, slice_index_y );

          sprintf( file, "%sP_Var_2D_XZ_j=%d.out", Dir );
          ((CPPEqn*) CPPe)->output_PVari( file, flag, slice_index_y );
      }

      cout << "Enter slice index in z-dir:... ";
      cin >> slice_index_z;
      if( slice_index_z >= 0 && slice_index_z < gridPtr->getNz() ) {
          flag = 7;
          sprintf( file, "%sP0_2D_XY_k=%d.out", Dir, slice_index_z );
          P0e->output( file, flag, *contPtr, slice_index_z );
          
          sprintf( file, "%sP2_2D_XY_k=%d.out", Dir, slice_index_z );
          P2e->output( file, flag, *contPtr, slice_index_z );
      
          sprintf( file, "%sCYP_2D_XY_k=%d.out", Dir, slice_index_z );
          CYPe->output( file, flag, *contPtr, slice_index_z );
      
          sprintf( file, "%sCPP_2D_XY_k=%d.out", Dir, slice_index_z );
          CPPe->output( file, flag, *contPtr, slice_index_z );

          sprintf( file, "%sP_Var_2D_XY_k=%d.out", Dir, slice_index_z );
          ((CPPEqn*) CPPe)->output_PVari( file, flag, slice_index_z );
      }
      */
  }
}

void Flow::writeToDomain(Domain * domainPtr) {
   domainPtr->setVelocity( veloPtr->getV0(0), veloPtr->getV0(1) );
   domainPtr->addVeloCova( veloPtr->getCv1v1(), veloPtr->getCv2v2(), veloPtr->getCv1v2() );
   domainPtr->calcDerivative();
}

