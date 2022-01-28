/*
 * File: Common.h
 * ------------------------------------------------------------------
 * This is a Utility file which will be used by everyone
 * 1) This file define all of the Constants and Enums.
 * 2) This file include all of the utility methods
 *
 * ------------------------------------------------------------------
 * Created    07/06/00          hcao     
 */

#ifndef _COMMON_H
#define _COMMON_H
#include <fstream>
#include <math.h>

using namespace std;

/* --- macro for machine difference ------------------- */
// --- Fortran call difference and bool type difference ---------
#ifdef SGI                  // SGI machine (SGI is defined by makefile)
#include <generic.h>
#define F77NAME(x) name2(x,_)

#else                    // IBM(SP2) machine
#define F77NAME(x) x

//#define bool int
//#define true 1
//#define false 0

#endif

/* --- first define the contants ---------------------- */
// the unit conversion factors (_ for "to"; P for "per"; 2 for "square";...)
const double ft_m    = 0.3048;
const double md_m2   = 1.0E-15 / 1.013;
const double day_sec = 24 * 3600;
const double stbPday_m3Psec = 5.615*pow(ft_m, 3)/day_sec;
const double psia_pa  = 1.013E+5/14.65;
const double cp_pasec = 0.001;
const double lbmPft3_kgPm3 = 0.454/pow(ft_m, 3);

// the universal contants 
const double PI = atan(1.0)*4.0;   


/* --- then define the enums for more physical meaning ------ */ 
enum Unit {
  FIELD  = 1,
  METRIC = 2
};

enum GridType {
  POINT_DIS = 1,
  CELL_CEN  = 2
};

enum WConType {        // well control type
  RATE = 1,
  PRES = 2
};

enum BType {          // boundary type
  CONST_RATE = 1,
  CONST_PRES = 2
};  

enum SolverType {
  FULLLU          = 1,
  LAPACK_BAND     = 2,
  LAPACK_FULL     = 3,
  LAPACK_SYM_BAND = 4
  //CLUB_SOLVER = 5
};

/* --- finally the utility methods (helper functions) ----------- */
// --- read and ignore the whole line after a certain character (# or /) ---
// --- used for comment in the input file (#... or //.....)
extern void Junk(ifstream &is);

// --- used for read an double array (one value for all, or discrete values)
extern double flagReading(ifstream &is, int n, double *x);

// --- used to find maximum of absolute value of an array (Pipatl)
extern double MaxAbs(int asize, double* arr);

#endif
