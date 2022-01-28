/*
 * File: ClubSolver.cc
 * ----------------------------------
 * Implementation for ClubSolver class
 */
#include "ClubSolver.h"

ClubSolver::ClubSolver(const Eqn &e) : Solver(e) 
{
  cout <<"Club Solver Constructor!" << endl;
  read_club_param();
  lrperm = lrtemp = litemp = 0;
  liperm = 1000;
  initiliaze( e );
}

ClubSolver::~ClubSolver() {
  delete [] iperm;
  delete [] ta;
  //club
  delete [] mgirow;
  delete [] mgjcol;
  delete [] ivec;
  delete [] tg;
  delete [] th;
  delete [] g;
  delete [] h;
  delete [] iws;
  delete [] dw;
  delete [] wrhs;
  delete [] pv;
  delete [] xcen;
  delete [] xoff;
  delete [] xrsd;
  delete [] rgtnsc;
  delete [] txoff;
  delete [] itemp;
  delete [] rperm;
  delete [] rtemp;
}

void ClubSolver::read_club_param() {
  ifstream in_param;
  in_param.open("club_param.dat",ios::in);
  if( !in_param ) {
      cerr << "Can not open club_param.dat file." << endl;
      exit(8);
  }
  
  in_param >> neqslv >> ipcont >> izcont >> iflag >> iorit  >> irsflg >> isp;
  in_param >> ndir;
  in_param >> nprm   >> nflowk >> north  >>mwells>> mlayx >> mcompr;
  //club
  in_param >> inorms >> jpcoar >> id4p   >> nit   >> iordin >> iludc;
  in_param >> restol;
  
  cout<<"Parameters of sqclub or club " << endl;
  cout << "neqslv = "<< neqslv << endl;
  cout << "ipcont = "<< ipcont << endl;
  cout << "izcont = "<< izcont << endl;
  cout << "iflag  = "<< iflag  << endl;
  cout << "iorit  = "<< iorit  << endl;
  cout << "irsflg = "<< irsflg << endl;
  cout << "isp    = "<< isp    << endl;
  cout << "nx     = "<< nx     << endl;
  cout << "ny     = "<< ny     << endl;
  cout << "nz     = "<< nz     << endl;
  cout << "ndir   = "<< ndir   << endl;
  cout << "nprm   = "<< nprm   << endl;
  cout << "nflowk = "<< nflowk << endl;
  cout << "north  = "<< north  << endl;
  cout << "mwells = "<< mwells << endl;
  cout << "mlayx  = "<< mlayx  << endl;
  cout << "mcompr = "<< mcompr << endl;
  cout<<"Parameters of club only" << endl;
  cout << "inorms = "<< inorms << endl;
  cout << "jpcoar = "<< jpcoar << endl;
  cout << "id4p   = "<< id4p   << endl;
  cout << "nit    = "<< nit    << endl;
  cout << "iordin = "<< iordin << endl;
  cout << "iludc  = "<< iludc  << endl;
  cout << "restol = "<< restol << endl;
  cout<<"Parameters of derived from given ones" << endl;          
  nflow =  ( nx + (nx+1)%2 )
          *( ny + (ny+1)%2 )
          *( nz + (nz+1)%2 ) + 1 ;
  ntrans= nx * ny * nz;
  nprmq = nprm * nprm;
  if(iordin == 0) nprrp = - nprm;
  else            nprrp =   nprm;
  nflt = nflowk;
  if(nflowk == 0) nflowk = 1;
  
  cout << "nflow  = "<< nflow  << endl;
  cout << "ntrans = "<< ntrans << endl;
  cout << "nprmq  = "<< nprmq  << endl;
  cout << "nprrp  = "<< nprrp  << endl;
  cout << "nflt   = "<< nflt   << endl;
  cout << "nflowk = "<< nflowk << endl;

  in_param.close();
}

void ClubSolver::initiliaze(const Eqn &e) {
  sqclub_init();
   setup_ta( e );
   calcLiperm();
   //exit(0);
    club_init();
}

void ClubSolver::sqclub_init() {
  ta = new double [ nflow * ndir ];
  int i;
  for(i = 0; i < nflow * ndir; ++i) ta[i] = 0.;
}

void ClubSolver::setup_ta(const Eqn &e) {
   int i, j, l;
   if( nx > 1) {                                 // assign off-D X+ and X-
       for(l = 0; l < nprmq; ++l) {
           for(i = 0; i < ntrans; ++i) {
               ta[i + l * nflow + 0 * nflow * nprmq] = e.offDAXp[i];
           }
       }
   }
   if( ny > 1) {                                   // assign off-D Y+ and Y-
       for(l = 0; l < nprmq; ++l) {
           for(i = 0; i < ntrans; ++i) {
               ta[i + l * nflow + 1 * nflow * nprmq] = e.offDAYp[i];
           }
       }
   }
  
   if( nz > 1) { 
       int nxy = nx*ny;  // assign off-D Z+ and Z-
       for(l = 0; l < nprmq; ++l) {
           for(i = 0; i < ntrans; ++i) {
               ta[i + l * nflow + 2 * nflow * nprmq] = e.offDAZp[i];
           }
       }
   }
}

void ClubSolver::calcLiperm() {
   int itry;
   int iperm_size;
   for(itry = 0; itry < 10; ++itry) {
       cout << "itry = " << itry <<endl;
       iperm = new int[ liperm ];
       cout << "iperm = new[]\n";
       iperm_size = liperm;
       SqClub_cc();

       cout << " lrperm = " << lrperm ; 
       cout << " liperm = " << liperm ;
       cout << " lrtemp = " << lrtemp ;
       cout << " litemp = " << litemp << endl;
       if ( liperm > iperm_size ) {
          cout << " iperm_size = " << iperm_size;
          cout << " lrperm = " << lrperm ; 
          cout << " liperm = " << liperm ;
          cout << " lrtemp = " << lrtemp ;
          cout << " litemp = " << litemp << endl;
          delete [] iperm;   
          cout << "delete[] iperm\n";  
       } else {
          break;
       } 
   }
}

void ClubSolver::club_init() {
  int i;

  mgirow = new int [ nflowk ];
  mgjcol = new int [ nflowk ];
  for(i = 0; i < nflowk; ++i) mgirow[i] = mgjcol[i] = 0;
  ivec   = new int    [ mcompr ];
  tg     = new double [ mcompr ];
  th     = new double [ mcompr ];
  for(i = 0; i < mcompr; ++i) {
     tg[i]   =  th[i]   =  0.;
     ivec[i] = 0;
  }
  g      = new double [ nprm * mcompr ];
  h      = new double [ nprm * mcompr ];
  for(i = 0; i < nprm * mcompr; ++i) g[i] = h[i] = 0.;
  dw     = new double [ mwells ];
  wrhs   = new double [ mwells ];
  for(i = 0; i < mwells; ++i) dw[i] = wrhs[i] = 0.;
  iws    = new int    [ mwells + 1 ];
  for(i = 0; i < mwells + 1; ++i) iws[i]  = 0;
  pv     = new double [ ntrans ];
  for(i = 0; i < ntrans ; ++i) pv[i] = 0.;
  xcen   = new double [ nflow * nprmq ];
  for(i = 0; i < nflow * nprmq; ++i) xcen[i] = 0.;
  xoff   = new double [ nflow  * nprmq * 2 * ndir ];
  for(i = 0; i < nflow  * nprmq * 2 * ndir; ++i) xoff[i] = 0.;
  xrsd   = new double [ nflow  * nprm ];
  for(i = 0; i < nflow  * nprm; ++i) xrsd[i] = 0.;
  rgtnsc = new double [ nflowk * nprmq ];
  for(i = 0; i < nflowk * nprmq; ++i) rgtnsc[i] = 0.;
  txoff  = new double [ nflow * 2 * ndir ];
  for(i = 0; i < nflow * 2 * ndir; ++i) txoff[i] = 0.;
  itemp  = new int    [ litemp ]; 
  for(i = 0; i < litemp; ++i) itemp[i] = 0;
  rperm  = new double [ lrperm ];
  for(i = 0; i < lrperm; ++i) rperm[i] = 0.;
  rtemp  = new double [ lrtemp ];
  for(i = 0; i < lrtemp; ++i) rtemp[i] = 0.;
}

void ClubSolver::setupA(const Eqn &e) 
{
   int i, j, k, l, n;

   // all are active cells
   for(i = 0; i < ntrans; ++i) pv[i] = 1.0;
   
   //SET UP COEFFICIENT MATRIX -- diagonal
   for(l = 0; l < nprmq; ++l) {
      for(i = 0; i < ntrans; ++i) {
          xcen[i + l * nflow ] = e.diagA[i + l * nflow ];
      }
   }
   int nxy = nx * ny;
   if( nx > 1) { 
       for(l = 0; l < nprmq; ++l) {   
           for(k = 0; k < nz; ++k) {
               for(j = 0; j < ny; ++j) {
		   n = j * nx + k * nxy;
		   xoff[n + l * nflow + 1 * nflow * nprmq] = e.offDAXp[n    ];
	           for(i = 1; i < nx - 1; ++i) {
		       n = i + j * nx + k * nxy;
		       xoff[n + l * nflow + 0 * nflow * nprmq] = e.offDAXn[n - 1];
		       xoff[n + l * nflow + 1 * nflow * nprmq] = e.offDAXp[n    ];
	           }
		   n = nx - 1 +j * nx + k * nxy;
		   xoff[n + l * nflow + 0 * nflow * nprmq] = e.offDAXn[n - 1];
	       }
           }
       }
   }
   if( ny > 1) {
       for(l = 0; l < nprmq; ++l) {
           for(k = 0; k < nz; ++k) {
	       for(i = 0; i < nx; ++i) {
		   n = i + k * nxy;
		   xoff[n + l * nflow + 3 * nflow * nprmq] = e.offDAYp[n    ];
	       }
               for(j = 1; j < ny - 1; ++j) {
	           for(i = 0; i < nx; ++i) {
		       n = i + j * nx + k * nxy;
		       xoff[n + l * nflow + 2 * nflow * nprmq] = e.offDAYn[n - nx];
		       xoff[n + l * nflow + 3 * nflow * nprmq] = e.offDAYp[n    ];
	           }
	       }
	       for(i = 0; i < nx; ++i) {
		   n = i - nx + (k + 1) * nxy;
		   xoff[n + l * nflow + 2 * nflow * nprmq] = e.offDAYn[n - nx];
	       }
           }
       }
   }
   if( nz > 1) { 
       for(l = 0; l < nprmq; ++l) {	   
           for(j = 0; j < ny; ++j) {
	       for(i = 0; i < nx; ++i) {
	           n = i + j * nx;
		   xoff[n + l * nflow + 5 * nflow * nprmq] = e.offDAZp[n    ];
	       }
	   }
           for(k = 1; k < nz - 1; ++k) {
               for(j = 0; j < ny; ++j) {
	           for(i = 0; i < nx; ++i) {
		       n = i + j * nx + k * nxy;
		       xoff[n + l * nflow + 4 * nflow * nprmq] = e.offDAZn[n - nxy];
		       xoff[n + l * nflow + 5 * nflow * nprmq] = e.offDAZp[n    ];
	           }
	       }
           }
	   k = nz - 1;
           for(j = 0; j < ny; ++j) {
	       for(i = 0; i < nx; ++i) {
	           n = i + j * nx + k * nxy;
		   xoff[n + l * nflow + 4 * nflow * nprmq] = e.offDAZn[n - nxy];
	       }
	   }
       }
   }
   
   
   // well part
   /*
   iws[0] = 1;
   for(i = 0; i < mwells; ++i) {
       dw[i] = 1.;
       iws[i + 1] = iws[ i ] + 1 ;
       wrhs[i] = 0.;
   }
   for(j = 0; j < mcompr; ++j) {
         tg[j] = 0.;
         th[j] = 0.;
         ivec[j] = iw[ j ] + jw[ j ] * nx + 1;
   }
   */
   // debug part
   /*
   ofstream os("club.debug", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(6);
   for(i = 0; i < ntrans; ++i) {
       os << xcen[i]<<' '
          << xoff[i + 0 * nflow] <<' '
	  << xoff[i + 1 * nflow] <<' '
	  << xoff[i + 2 * nflow] <<' '
	  << xoff[i + 3 * nflow] <<' '
	  << xrsd[i] << endl;
   }
   os.close();
   */
}

void ClubSolver::solve(bool doLU, const Eqn &e, int nRHS, double *RHS)
{
  int i, j, l;
  double startT, endT;
  //cout << nx <<' '<< ny <<' '<< nz <<' '<<ntrans<< endl;
  for(l = 0; l < nRHS; ++l) {
      for(j = 0; j < nprm; ++j) {
          for(i = 0; i < ntrans; ++i) {
              xrsd[ i + j * nflow ] = RHS[l * ntrans + i];
          }
      }
   
      if( doLU ) { 
          setupA( e );  
          //startT = clock();
          Club_cc();
          //endT   = clock();
	  //cout << "l = " << l<<" ntns = " << ntns << endl;
	  if(ntns >= 100) {
	     cout << "l = " << l<<" ntns = " << ntns << endl;
	     exit(0);
	  }
      } else {
          iflag = - 1;
	  Club_cc();
	  if(ntns >= 100) {
	     cout << "l = " << l<<" ntns = " << ntns << endl;
	     exit(0);
	  }
	  //cout << "l = " << l<<" ntns = " << ntns << endl;
      }

      for(i = 0; i < ntrans; ++i) {
          RHS[l * ntrans + i] =   xrsd[ i ];
      }
  }
}

/* --- Private Methods ----------------------------- */
// --- C++ requirement to call Fortran subroutine ----
extern "C"
{
  void sqclub_(int* neqslv, int* ipcont, int* izcont, int* iflag , int* iorit ,
               int* irsflg, int* isp   , int* nx    , int* ny    , int* nz    ,
               int* ndir  , int* nprm  , int* nflow , double* ta , int* nflowk,
               int* north , int* mwells, int* mlayx , int* mcompr, int* iperm ,
               int* lrperm, int* liperm, int* lrtemp, int* litemp );

  void   club_(int* inorms, int* jpcoar   , int* id4p  , int* ipcont, int* izcont, 
               int* north , double* restol, int* nit   , int* iordin, int* nprm  , 
               int* nprmq , int*    nprrp , int* iludc , int* iorit , int* ntns  , 
               int* irsflg, int*    ndir  , int* mcompr, int* mlayx , int* neqslv, 
               int* isp   , int*    nx    , int* ny    , int* nz    , int* ndir  , 
               int* ntrans, int*    nflow , int* nflowk, int* nflt  , int* mgirow, 
               int* mgjcol, int*    mwells, double*  g , double* h  , double* dw , 
               double* wrhs, double* pv    , double*xcen, double* xoff, double* xrsd, 
               double* rgtnsc, int* iws    , int * ivec , double* txoff,double* tg,
               double* th , int* iperm , double* rperm, int* itemp,double* rtemp);
};

// --- Wrappers to SQClub Fortran subroutine ------
void ClubSolver::SqClub_cc()
{
  sqclub_(&neqslv, &ipcont, &izcont, &iflag , &iorit ,
          &irsflg, &isp   , &nx    , &ny    , &nz    ,
          &ndir  , &nprm  , &nflow ,  ta    , &nflowk,
          &north , &mwells, &mlayx , &mcompr,  iperm ,
          &lrperm, &liperm, &lrtemp, &litemp);
}

// --- Wrappers to Club Solver Fortran subroutine  -------------
void ClubSolver::Club_cc()
{ 
  /*
  cout << inorms << ' '<<jpcoar<<' '<<id4p <<' '<<ipcont <<' '<<izcont<<endl;
  cout << north  << ' '<<restol<<' '<<nit  <<' '<<iordin <<' '<<nprm  <<endl;
  cout << nprmq  << ' '<<nprrp <<' '<<iludc<<' '<<iorit  <<' '<< ntns <<endl;
  cout << irsflg << ' '<<ndir  <<' '<<mcompr<<' '<<mlayx <<' '<< neqslv<< endl;
  cout << isp    << ' '<<nx    <<' '<<ny    <<' '<<nz    <<' '<< ndir  << endl;
  cout << ntrans << ' '<<nflow <<' '<<nflowk<<' '<<nflt  <<' '<< mwells<<endl;
  */
  club_(&inorms, &jpcoar, &id4p  , &ipcont, &izcont,
        &north , &restol, &nit   , &iordin, &nprm  ,
        &nprmq , &nprrp , &iludc , &iorit , &ntns  ,
        &irsflg, &ndir  , &mcompr, &mlayx , &neqslv,
        &isp   , &nx    , &ny    , &nz    , &ndir  ,
        &ntrans, &nflow , &nflowk, &nflt  ,  mgirow,
         mgjcol, &mwells,  g     ,  h     ,  dw    ,
         wrhs  ,  pv    ,  xcen  ,  xoff  ,  xrsd  ,
         rgtnsc, iws    ,  ivec  ,  txoff ,  tg    ,
         th    , iperm  ,  rperm ,  itemp ,  rtemp );

  //cout << "club_cc() is used, the iteration number = " << ntns << endl;
  /*
  cout << "after\n";
  cout << inorms << ' '<<jpcoar<<' '<<id4p <<' '<<ipcont <<' '<<izcont<<endl;
  cout << north  << ' '<<restol<<' '<<nit  <<' '<<iordin <<' '<<nprm  <<endl;
  cout << nprmq  << ' '<<nprrp <<' '<<iludc<<' '<<iorit  <<' '<< ntns <<endl;
  cout << irsflg << ' '<<ndir  <<' '<<mcompr<<' '<<mlayx <<' '<< neqslv<< endl;
  cout << isp    << ' '<<nx    <<' '<<ny    <<' '<<nz    <<' '<< ndir  << endl;
  cout << ntrans << ' '<<nflow <<' '<<nflowk<<' '<<nflt  <<' '<< mwells<<endl;
  */
}

void ClubSolver::print() {
   ofstream os("club_driver1.out", ios::out);
   os << setiosflags(ios::fixed | ios::showpoint)
      << setprecision(6);
   int i, j;
   for(j = 0; j < ny; ++j) {
      for(i = 0; i < nx; ++i) {
         os << xrsd[ i + j * nx ] << endl;
      }
      os << endl;
   }
   os.close();
}

