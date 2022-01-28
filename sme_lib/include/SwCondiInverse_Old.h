/*
 * File: SwCondiInverse.h
 * This file define the SwCondiInverse class
 * Auther: Pipat Likanapaisal
 */

#ifndef _SwCondiInverse_h
#define _SwCondiInverse_h

#define NO_INDEX -1

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <iomanip>
#include <math.h>

#include "Grid.h"
#include "Point.h"
#include "Domain.h" //This is necessary for backtracking only(for now).
#include "ParticleTrack.h" //This is necessary for backtracking only(for now).
#include "Perm.h"
#include "Control.h"
#include "Solver.h"
#include "Well.h"
#include "Flow.h"
#include "BLSolution.h"
#include "Saturation.h"
#include "KrigWeigh.h"
//#include "quadrature.cc"

void GaussianQuadrature(int qPoint, double* &rloc, double* &weight);

/* Class: SwCondiInverse
 * ---------------------
 * This class is used to store transport conditioning data.
 */

class SwCondiInverse
{
	public:
		SwCondiInverse(Domain* domain, Flow* flow);
		~SwCondiInverse();

		void SaturationConditioning(char *directory, ifstream &is_flow);
		void RefreshPressureGradient();
		double CalcCovYVx(int i1, int j1, int k1, int i2, int j2, int k2, int yInd);
		double CalcCovYVy(int i1, int j1, int k1, int i2, int j2, int k2, int yInd);
		//void testBackTracking();
		
		//For Production to use to calculate water production rate
		void CollectDataS(int lInd);
		void DeallocateDataS();		
		double CalcCovYTau(int lInd);
		double CalcCovYLnTau(double tauMean, double covYTau);
		void CalcCov_YVxi_YVeta_YEta(int lInd, int yInd, int pInd);
		void BTrackPtrIs(ParticleTrack* btp){bTrackPtr = btp;}

	private:
		struct Coord;
		//struct SwCondData;

		void ReadData();
		Coord* FindCoordinate(Coord* &coordPtr, Coord &coordVal);
		int CompareCoordinate(Coord coord1, Coord coord2);
		void AddCoordinate(Coord* &coordPtr, Coord coordVal);
		void SequentialConditioning(char *directory, ifstream &is_flow);
		void GlobalConditioning(char *directory, ifstream &is_flow);
		void ResetCoordSetIndex(Coord* coordPtr);
		int CountMeasurements(int seqLoc);
		//void LaunchParticles(int seqLoc, int numSw);
		void LaunchParticles(Coord** condPtr, int numCond, Coord** obvPtr, int numObv, int iter);
		void DisplayParticleInfo(int npoint, ParticleTrack* pTrackPtr);
		void CalcSaturationMoments(Coord** cPtrPtr, int pLoc, int numP, ofstream &pof, int iter, int condPointFlag);
		void DeallocateIter();
		void CalcCovSYVector(int iter, int numSw, int pLoc);
		void FindEdge(int lInd, int pInd);
		double CalcCovSYQuad(int sInd, int yInd, int pLoc); //Does similarly to CalcCovSY, but with quadrature integration
		void CalcLambdaVector(int numSw);
		void UpdatePermeabilityMoments(int seqLoc, int numSw, int iter);
		////CalcCovSY() is replaced by CalcCovSYQuad which has quadrature integration
		//double CalcCovSY(int sInd, int yInd);
		////Since the relation between CYTau and CYLnTau has been derived, the two functions below are not necessary
		//double RecCalcCovYLnTau(double yMean, double yVar, double tauMean, double tauVar, double covYTau,int ns, double eps, double covYTauMin, double covYLnTauMin, double covYTauMax, double covYLnTauMax);
		//double CalcCovYTauByIntegration(double yMean, double yVar, double tauMean, double tauVar, double covYLnTau, int ns);

		struct Coord
		{
			int i_,j_,k_;
			double x_,y_,z_;
			int index;
			Coord* coordPtr;
		};
		/*struct SwCondData
		{
			double sw;
			Coord* coordPtr;
		};*/

		Coord* coordSet;
		Coord** condSet;
		Coord** obvSet;

		int condScheme; // 0 = Sequential, 1 = Global
		double dFactor; // damping factor (multiplied to lambda)
		int nCond, nObv;
		double *condTime, *obvTime;
		double *swMeas;
		int *condSwStrl, *obvSwStrl;
		int nCoord;
		int totalCondiTime;
		int maxIter;
		double etol,vtol; //error tolerance, variance tolerance

		bool debug;
		double* satCovEP;
		double* satVarEP;
		double* satAvgEP;
		double* satMeas;
		double* satErr;
		double* covSY;
		double* covYVxiU;double* covYVxiD; // U = upstream, D = downstream
		double* covYVetaU;double* covYVetaD;
		double* covYEta;
		double* lambda;
		double  Swc, Sor, viscosR;
		int nS; //number of discretized sections for integration
		double clockStart, clockEnd;
		int readDataFlag;

		double* ds;
		double* dvdeta;
		double* costh;
		double* sinth;
		int* edge;
		double *x_, *y_, *z_;
		int *i_, *j_, *k_;
		double *xlb_, *ylb_, *zlb_; // Left bottom
		double *xrt_, *yrt_, *zrt_; // Right top
		const double* dpdx;
		const double* dpdy;

		Domain* domainPtr;
		Flow* flowPtr;
		Grid* gridPtr;
		Perm* permPtr;
		P0Eqn* P0e;
		CYPEqn* CYPe;

		Point* pointPtr;
		ParticleTrack* bTrackPtr;
		BLSolution* blSolnPtr;
		Saturation* condSatPtr;
		Saturation* obvSatPtr;
		KrigWeigh* krigPtr;
};
#endif

