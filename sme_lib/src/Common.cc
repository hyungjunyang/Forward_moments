// --- implementation to utility functions (C Style, no class) ---
#include "Common.h"
#include <fstream>

// --- read and ignore the whole line after a certain character (# or /) ---
// --- used for comment in the input file (#... or //.....)
void Junk(ifstream &is) {
  char sym;
  while(true) {
    is >> sym;
    //    if(sym!='#' && sym!='/') break;
    if(sym!='#') break;
    is.ignore(128, '\n');  // ignore the comment line
  }
  is.putback(sym);    // put back the non-comment character
}

// --- used for read an double array (one value for all, or discrete values)
double flagReading(ifstream &is, int n, double *x) {
  int i, flag; double tmp;
  Junk(is); is >> flag; 
  Junk(is);
  if(flag == 1) {    // one value for all
    is >> tmp; 
    for( i = 0; i < n; i++) x[i] = tmp;
  } 
  else {         // discrete values 
    for( i = 0; i < n; i++) is >> x[i];
  }
  return tmp;  // return the single(one) value
}

// --- used to find maximum of absolute value of an array (Pipatl)
extern double MaxAbs(int asize, double* arr)
{
	double maxValue = fabs(arr[0]);
	for(int ii=1;ii<asize;ii++)
	{
		if(fabs(arr[ii])>maxValue){maxValue = fabs(arr[ii]);}
	}
	return maxValue;
}

