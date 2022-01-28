/*
 *
 * --------------------------------------------------------------------
 * Created    08/05/03          liyl     
 */


#ifndef _CLUBSOLVER_H
#define _CLUBSOLVER_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "Solver.h"
#include "Eqn.h" 

using namespace std;

class ClubSolver : public Solver{

 public:
   ClubSolver(const Eqn &e);
   virtual ~ClubSolver();
 
   // --- Multi RHS direct solve -----------------
   virtual void solve(bool doLU, const Eqn &e, int nRHS, double *RHS);
   void print(); 
 protected:
   virtual void setupA(const Eqn &e); 

 private:
  /* sqclub variables*/
  int neqslv, ipcont, izcont, iflag , iorit;
  int irsflg, isp   ;
  int ndir  , nprm  , nflow , nflowk;
  int north , mwells, mlayx , mcompr; 
  int lrperm, liperm, lrtemp, litemp;
  int   * iperm;
  double* ta;

  /** club variables */
  int inorms, jpcoar, id4p  , nit   , iordin,  nprmq;
  int nprrp , iludc , ntns  , ntrans, nflt;
  double restol;
  
  int *mgirow, *mgjcol;
  int *itemp , *ivec, *iws;
  double *g,  *h, *dw, *wrhs, *pv;
  double *xcen, *xoff, *xrsd, *rgtnsc;
  double *txoff, *tg,  *th, *rperm, *rtemp;

  void read_club_param();
  void initiliaze(const Eqn &e);
  void sqclub_init();
  void club_init();
  void calcLiperm();
  void setup_ta(const Eqn &e);
  	  
  // --- Wrappers to Club PGF90 subroutines (sqclub, club) ------
  void SqClub_cc();
  void Club_cc(); 
};

#endif

/*
 * File: Club.h
 * --(SQClub)-----------------------------------------------------------------
     NEQSLV: NUMBER OF EQUATIONS PER GRID BLOCK
     IPCONT: FLAG OF PRESSURE CONSTRAINT FOR MATRICES with NEQSLV > 1
             =1, PRESSURE CONSTRAINT ON , RECOMMENDED FOR  NEQSLV > 1 
               (Not Available for MNF -- Modified Nested Factorization)
             =0, PRESSURE CONSTRAINT OFF, RECOMMENDED FOR  NEQSLV = 1
     IZCONT: FLAG OF Z-LINE CONSTRAINT,
             =1, Z-LINE CONSTRAINT ON,  RECOMMENDED FOR 3-D cases
             =0, Z-LINE CONSTRAINT OFF, RECOMMENDED FOR 2-D and 1-D
     IFLAG: FLAG FOR Storing and retrieving COMMON block variables in the club solver.
            = 1, Estimate memory requirements for club. 
                 The required array sizes are 
                 LIPERM, LRPERM, LITEMP, and LRTEMP.
            = 0, COPY COMMON TO IPERM, INITIALLY
                 IF IORDIN IS TO BE SET TO 0 IN NEXT CALL TO CLUB,
                 SET IFLAG = 0
            =-1, COPY IPERM TO COMMON, OTHERWISE
     IORIT: DIRECTION FOR ORDERING, RECOMMENDED VALUE = 0
         0      AUTOMATIC ORIENTATION (IT CHECKS ALL 3 DIRECTIONS FOR
                THE LEAST NUMBER OF BLOCKS)
         1      ZYX
         2      ZXY
         3      YXZ
         4      YZX
         5      XYZ
         6      XZY
        -1      DO NOT USE AXYZ, THIS IS FOR MEMORY ESTIMATION ONLY
     IRSFLG: FLAG FOR PRECONDITIONER
         0      DIAGONAL SCALING FOR IMPLICIT SOLVE, RS/ILU(0) FOR IMPES
                0 IS RECOMMENDED
         1      TRI-DIAGONAL SCALING FOR IMPLICIT SOLVE, RS/ILU(0) FOR IMPES
         3      MNF for Implicit, RS/ILU(0) for IMPES
         5      DIAGONAL SCALING for Implicit, MNF for IMPES
         6      TRI-DIAGONAL SCALING for Implicit, MNF for IMPES
         8      MNF for Implicit, MNF for IMPES
     ISP   : FLAG FOR SPECIAL MATRIX STRUCTURE
         0      7-POINT, SINGLE PERMEABILITY
         1     11-POINT, SINGLE PERMEABILITY
         2     7-POINT, DUAL-PERMEABILITY
         3     11-POINT, DUAL-PERMEABILITY
     NX,NY,NZ: DIMENSION OF 3-D MODEL
     MDIR: NUMBER OF DIMENSIONS (3 for 3-D)
     NPRM: NUMBER OF BALANCE EQUATIONS, USUALLY = NEQSLV
     NFLOW: FIRST DIMENSION OF XCEN, XOFF, XRSD, AXYZ
            SHOULD BE >= NX*NY*NZ
     TA(NFLOW,*): TRANSMISSIBILITY ARRAY
     NFLOWK: FIRST DIMENSION OF RGTNSC
             SHOULD BE >= THE NUMBER OF IRREGULAR CONNECTIONS
     NORTH: NUMBER OF ORTHOGONALIZATION, RECOMMENDED VALUE = 10
     MWELLS: NUMBER OF WELLS, SHOULD BE >= NWELLS
     MLAYX : MAXIMUM NUMBER OF COMPLETIONS PER WELL
     MCOMPR: TOTAL NUMBER OF COMPLETIONS FOR THE ENTIRE MODEL
     IPERM(): PERMANENT ARRAYS FOR INTEGER
            NEEDS AT LEAST 1000, IF LIPERM IS STILL UNKNOWN
     LRPERM: OUTPUT, DIMENSION FOR REAL PERMANENT ARRAYS
     LIPERM: OUTPUT, DIMENSION FOR INTEGER PERMANENT ARRAYS
     LRTEMP: OUTPUT, DIMENSION FOR REAL TEMPORARY ARRAYS
     LITEMP: OUTPUT, DIMENSION FOR INTEGER TEMPORARY ARRAYS
 * --(Club)------------------------------------------------------------
 *   INORMS : 0, UNNORMALIZED COEFFICIENTS ARE ENTERED
 *            1,   NORMALIZED COEFFICIENTS ARE ENTERED
 *   JPCOAR : 0, SOLVE FULLY IMPLICIT MATRIX
 *            1, SOLVE PRESSURE CONSTRAINT MATRIX ONLY (NO FULLY IMPLICIT MATRIX)
 *   ID4P   : POSITION OF THE PRESSURE UNKONW IN NEQSLV
 *   IPCONT, IZCONT, NORTH, SEE SQCLUB
 *   RESTOL: TOLERANCE, RECOMMENDED VALUE = 10**-4 FOR PRESSURE 
 *                                        = 10**-3 FOR MULTIPLE UNKNOWS
 *   NIT   : MAXIMUM NUMBER OF ITERATIONS, SHOULD BE LESS THAN 50
 *   IORDIN: FLAG FOR ORDERING MATRIX
 *           = 0, ORDERING MATRIX INDEX, 
                  SHOULD BE SET TO O EVERYTIME INCIDENT MATRIX IS CHANGED
             = 1, DO NOT REORDER MATRIX INDEX
     NPRM, SET SQCLUB
     NPRMQ : NPRM*NPRM
     NPRRP : = -NPRM, IF IORDIN = 0
               NPRM,  IF IORDIN = 1
     ILUDC : FLAG FOR ILU DECOMPOSITION
             = 1, DECOMPOSE, MOST OF THE TIME USE THIS VALUE
             = 0, DO NOT DECOMPOSE
     IORIT : SEE SQCLUB
     IRSFLG, ISP,  MDIR, MCOMPR, MLAYX, NEQSOL(=NEQSLV), 
     NX, NY, NZ, NDIR(=MDIR): SEE SQCLUB
     NTRANS: NX*NY*NZ
     NFLOW, NFLOWK: SEE SQCLUB
     NFLT  : NUMBER OF IRREGULAR CONNECTIONS <= NFLOWK
     MGIROW, MGJCOL: IRREGULAR CONNECTION INDEXES
     NWELLS: NUMBER OF WELLS
     G, H, DW, WRHS, XCEN, XOFF, XRSD, RGTNSC, 
               IWS, IVEC, FOLLOW BLITZ NOTATIONS
     XRSD  : RIGHT HAND SIDE RESIDUAL, ON OUTPUT, IT STORE SOLUTION
     TXOFF : DIMENSION TO NFLOW*NDIR*2, TEMPORARY ARRAY, NEEDED ONLY
             IF XOFF STORED AS XOFF(NFLOW,NPRM,*) IN THE 1-EQUATION CASE
     TG    : DIMENSION TO MCOMPR, TEMPORARY ARRAY, NEEDED ONLY
             IF G STORED AS G(NPRM,*) IN THE 1-EQUATION CASE
     TH    : DIMENSION TO MCOMPR, TEMPORARY ARRAY
             IF G STORED AS H(NPRM,*) IN THE 1-EQUATION CASE
     IPERM : SIZE OF LIPERM, PERMANENT ARRAY
     RPERM : SIZE OF LRPERM, PERMANNET ARRAY
     ITEMP : SIZE OF LITEMP, TEMPORARY ARRAY
     RTEMP : SIZE OF LRTEMP, TEMPORARY ARRAY

 *   OUTPUt
 *   NTNS            : Total number of iterations
 *   XRSD(NFLOW,NPRM): Solution vector
 *   WRHS(NWELLS)    : Bottom hole pressures solution vector.
 *   Iperm(Liperm)   : This array stores Common block variables and ILU
 *                     information.
 *   Rperm(Lrperm)   : This array stores ILU information.
 * --------------------------------------------------------------------
 * Created    07/29/03          liyl     
 * */
