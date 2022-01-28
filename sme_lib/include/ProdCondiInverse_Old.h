/*
 * File: ProdCondiInverse.h
 * This file define the ProdCondiInverse class
 * Auther: Pipat Likanapaisal
 */

#ifndef _ProdCondiInverse_h
#define _ProdCondiInverse_h

#define INT_SECTION_PER_SIDE	5;
#define INT_QUAD_POINT		3;

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <math.h>

#include "Common.h"
//#include "Grid.h"
//#include "Point.h"
#include "Domain.h"
//#include "ParticleTrack.h"
//#include "Perm.h"
//#include "Control.h"
//#include "Solver.h"
#include "Well.h"
#include "Flow.h"
#include "BLSolution.h"
#include "Production.h"
#include "KrigWeigh.h"
#include "SwCondiInverse.h"

void GaussianQuadrature(int qPoint, double* &rloc, double* &weight);

/* Class: ProdCondiInverse
 * -----------------------
 * This class is used to store production conditioning data.
 */

class ProdCondiInverse
{
	public:
		ProdCondiInverse(Domain* domain, Flow* flow);
		~ProdCondiInverse();
		void ProductionConditioning(char *directory, ifstream &is_flow);

		//void ProductionConditioning(char *directory, ifstream &is_flow);

	private:
		struct Coord;
		struct WCondData;

		void ReadData();
		void CalcProductionMoments(ofstream &pof, int iter, int condWellFlag);
		void DeallocateIter();
		void CalcCovQQMatrix();
		void CalcCovQYVector();
		void CalcLambdaVector();
		void UpdatePermeabilityMoments(int iter);

		/*struct Coord
		{
			int i_,j_,k_;
			double x_,y_,z_;
			int index;
			Coord* coordPtr;
		};*/
		struct WCondData
		{
			int wIndex;
			int calcQtFlag;
			int numQw;
			//Coord** wCoordPtr;
		};

		Domain* domainPtr;
		Flow* flowPtr;
		Well* wellsPtr;
		Grid* gridPtr;
		Perm* permPtr;
		KrigWeigh* krigPtr;
		Production* prodPtr;
		BLSolution* blSolnPtr;
		SwCondiInverse* swCondPtr;

		//Coord* coordSet;
		WCondData* condSet;
		WCondData* obvSet;
		double* condQtMeas;
		double** condQwTime;
		double** condQwMeas;
		double** obvQwTime;

		double dFactor; // damping factor (multiplied to lambda)
		int nWCond, nWObv;
		int nMeas;
		int nObv;
		int maxIter;
		double etol,vtol; //error tolerance, variance tolerance
		double clockStart, clockEnd;
		int readDataFlag;

		double* prodCov;
		double* prodVar;
		double* prodAvg;
		double* prodErr;
		double* covQY;
		double* lambda;

		double Swc, Sor, viscosR;
};
#endif

