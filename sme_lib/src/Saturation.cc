/*
 * File: Saturation.cc
 * ----------------------------------
 * Implementation for Saturation class
 */
#include "Saturation.h"

/* --- Public Methods ---------------------------------- */
// --- Constructors and destructor -----------
Saturation::Saturation(ParticleTrack *ptclTrackPtr_, 
		       double visR_, double Swc_, double Sor_, int nS_ )
 : ptclTrackPtr(ptclTrackPtr_), visR( visR_ ), Swc( Swc_ ), Sor( Sor_ ), nS ( nS_ )
{
  deltaS = 1.0 - Swc - Sor;

  // --- assign streamline data -----
  
  nLn = ptclTrackPtr->getNumStrl();
  nNd = 0; 
  for(int i = 0; i < nLn; ++i) {
      int len_tmp = ptclTrackPtr->getParticle(i)->getLength();
      if(len_tmp > nNd) nNd = len_tmp;
  }
 
  nPoint = nNd * nLn;
  
  // --- allocate saturation data ---
  initialize();
}

Saturation::~Saturation() {
  cout << "Saturation::~Saturation()" << endl;
  delete[] tau; 
  delete[] s;
  delete[] satDet; 
  delete[] satAvg; 
  delete[] satVar;
  delete[] tauMn;
  delete[] tauVar;
  delete[] x1C;
  delete[] etaMnC;
}

// --- initialize saturation data ----------------
void Saturation::initialize() {
  tau    = new double[nS];
  s      = new double[nS];
  tauMn  = new double[nPoint];
  tauVar = new double[nPoint];
  x1C    = new double[nPoint];
  etaMnC = new double[nPoint];	  
  for(int i = 0; i < nPoint; ++i) {
      tauMn[i] = tauVar[i] = 0; 
        x1C[i] = etaMnC[i] = 0;
  }
  for(int i = 0; i < nLn; ++i) {
      int len_tmp = ptclTrackPtr->getParticle(i)->getLength();
      for(int j = 0; j < len_tmp; ++j) {
          tauMn [j + i * nNd] = ptclTrackPtr->getParticle(i)->getTrvlTimeAvg(j);
          tauVar[j + i * nNd] = ptclTrackPtr->getParticle(i)->getTrvlTimeVar(j);
	      x1C[j + i * nNd] = ptclTrackPtr->getParticle(i)->getTrajectory(j).getX();
          etaMnC[j + i * nNd] = ptclTrackPtr->getParticle(i)->getTrajectory(j).getY();
      }
  }
}

// --- solver for saturation moments at time t---------
/*void Saturation::solve(double tt) {
  t = tt;
  //int iln, jnd, ij;
  double startT, endT;

  // --- calculate sStar, fpwStar ----------------
  sStar = (1-Sor)*(1-Sor)+visR*Swc*Swc-2*deltaS*(1-Sor);
  sStar /= visR+1;
  sStar = Swc+sqrt(Swc*Swc-sStar);  
  fpwStar = fpw(sStar);
  cout << "sStar=" << sStar << ",  fpwStar=" << fpwStar << endl;

  // --- calculate deterministic saturation ------
  startT = clock();
  calcSatDet(nPoint, tauMn, satDet);
  endT   = clock();
  cerr << "Sat Deter time: " <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;

  // --- calculate saturation mean and vari ------
  startT = clock();
  calcSatAvgVarQuad();
  endT   = clock();
  cerr << "Sat MnVari time: " <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;

}*/
void Saturation::solve(double tt)
{
	double* timeArray = new double[nLn];
	int* swStrl = new int[nLn];
	for(int ii=0;ii<nLn;ii++)
	{
		timeArray[ii] = tt;
		swStrl[ii] = ii;
	}
	solveMT(timeArray,nLn,swStrl,0);
	delete[] timeArray;
	delete[] swStrl;
}

// --- calculate saturation for multiple time
// timeArray = array of the multiple time
// tSize = size of timeArray
// swStrl = the number of streamline used to calculate saturation at timeArray
// endPointFlag: 0 = calculate saturation for all point on streamlines
//				 1 = calculate saturation for only the last point on streamlines
void Saturation::solveMT(double* timeArray, int tSize, int* swStrl, int endPointFlag)
{
	this->timeArray = timeArray;
	this->tSize = tSize;
	this->swStrl = swStrl;
	this->endPointFlag = endPointFlag;
	double startT, endT;

	int satSize = tSize*(((1-endPointFlag)*nNd)+(endPointFlag));
	satDet = new double[satSize]; 
	satAvg = new double[satSize];
	satVar = new double[satSize];

	// --- calculate sStar, fpwStar ----------------
	sStar = (1-Sor)*(1-Sor)+visR*Swc*Swc-2*deltaS*(1-Sor);
	sStar /= visR+1;
	sStar = Swc+sqrt(Swc*Swc-sStar);  
	fpwStar = fpw(sStar);
	cout << "sStar=" << sStar << ",  fpwStar=" << fpwStar << endl;

	// --- calculate deterministic saturation ------
	startT = clock();
	for(int istrl=0;istrl<tSize;istrl++)
	{
		int len_tmp = ptclTrackPtr->getParticle(swStrl[istrl])->getLength();
		int pointStart;
		if(endPointFlag==0){pointStart = 0;}
		else{pointStart = len_tmp-1;}
		t = timeArray[istrl];
		calcSatDet(len_tmp-pointStart,&tauMn[(swStrl[istrl]*nNd)+(endPointFlag*(len_tmp-1))],&satDet[istrl*(((1-endPointFlag)*nNd)+(endPointFlag))]);
	}
	endT   = clock();
	cerr << "Sat Deter time: " <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;

	// --- calculate saturation mean and vari ------
	startT = clock();
	calcSatAvgVarQuad();
	endT   = clock();
	cerr << "Sat MnVari time: " <<(endT-startT)/CLOCKS_PER_SEC << " sec" << endl;
}

// --- calculate deterministic solution -----------------
void Saturation::calcSatDet(int n, const double *tau1, double *satDet) {
  double eps = 1.0E-12;

  for(int ij = 0; ij < n; ij++) {
           if(tau1[ij] <=           0) satDet[ij] = 1.0 - Sor;               // inlet
      else if(tau1[ij] == fpwStar * t) satDet[ij] = sStar;                   // front
      else if(tau1[ij] >  fpwStar * t) satDet[ij] = Swc;                     // ahead front
      else                             satDet[ij] = root(tau1[ij]/t, eps);   // behind front
  }
}

void Saturation::calcSatAvgVarQuad()
{
	int ij, ijsat, iq;
	double sDRange = 3.5, intRange = 1;
	int qPoint = 3;
	double* rloc, * weight;
	GaussianQuadrature(qPoint,rloc,weight);

	// --- calculate saturation mean and variance ----------
	for(int istrl = 0; istrl < tSize; istrl++) {
		t = timeArray[istrl];
		int len_tmp = ptclTrackPtr->getParticle(swStrl[istrl])->getLength();
		int pointStart;
		if(endPointFlag==0){pointStart = 0;}
		else{pointStart = len_tmp-1;}
		for(int jnd = pointStart; jnd < len_tmp; jnd++) {
			ij = (swStrl[istrl]*nNd) + jnd;
			ijsat = (istrl*(((1-endPointFlag)*nNd)+(endPointFlag))) + ((1-endPointFlag)*jnd);
			if(tauMn[ij]==0) //tauMn == 0 ==> tauVar == 0
			{
				satAvg[ijsat] = satDet[ijsat];   // quite deterministic
				satVar[ijsat] = 0.0;
			}
			else
			{
				// do the integration
				double lntauMn,lntauVar,sVal;
				double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
				int lenn = 1;
				transform1(lenn,&tauMn[ij],&tauVar[ij],&lntauMn,&lntauVar);
				lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
				lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
				lntauFront = log(fpwStar*t);
				satAvg[ijsat] = satVar[ijsat] = 0.0;
				while(true)
				{
					rPoint = lPoint+(intRange*sqrt(lntauVar));
					if(rPoint>lntauMax){rPoint = lntauMax;}
					if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
					for(iq = 0; iq < qPoint; iq++) {
						double lntauVal = ((rPoint-lPoint)/2*rloc[iq])+((rPoint+lPoint)/2);
						double tauVal = exp(lntauVal);
						double sVal;
						calcSatDet(1,&tauVal,&sVal);
						double tmp = (rPoint-lPoint)/2*weight[iq]*sVal*pdfGauss(lntauVal,lntauMn,lntauVar);
						satAvg[ijsat]   += tmp;
						satVar[ijsat] += tmp*sVal;
					}
					if(rPoint==lntauMax){break;};
					lPoint = rPoint;
				}
				satVar[ijsat] -= (satAvg[ijsat]*satAvg[ijsat]);
			}
		}
	}
	delete[] rloc;
	delete[] weight;
}

double* Saturation::calcSatCovEndPoints()
{
	// Initialize
	double* satCovEP = new double[tSize*tSize];
	rhoLntau = ptclTrackPtr->getLogTrvlTimeRho();
	ofstream aof("Css_rhoLnTau.out", ios::out);
	aof << "icond1 icond2 Sw_Var rho_LnTau" << endl;
	for(int ii=0;ii<tSize;ii++)
	{
		for(int jj=ii;jj<tSize;jj++)
		{
			int iijj = ii*tSize + jj;
			int isnd = (ii*(((1-endPointFlag)*nNd)+(endPointFlag))) + ((1-endPointFlag)*(ptclTrackPtr->getParticle(swStrl[ii])->getLength()-1));
			if(ii==jj)
			{
				if(satVar[isnd]>0){satCovEP[iijj] = satVar[isnd];}
				else{satCovEP[iijj] = SMALL_NUM;} // This is to prevent covariance matrix become singular
				aof<<jj<<"\t"<<ii<<"\t"<<satVar[isnd]<<"\t"<<1<<endl;
			}
			else
			{
				int jsnd = (jj*(((1-endPointFlag)*nNd)+(endPointFlag))) + ((1-endPointFlag)*(ptclTrackPtr->getParticle(swStrl[jj])->getLength()-1));
				int itnd = (swStrl[ii]*nNd)+(ptclTrackPtr->getParticle(swStrl[ii])->getLength()-1);
				int jtnd = (swStrl[jj]*nNd)+(ptclTrackPtr->getParticle(swStrl[jj])->getLength()-1);
				int ijtnd = swStrl[ii]*nLn + swStrl[jj];
				int jjii = jj*tSize + ii;
				satCovEP[iijj] = calcSatCovQuad(timeArray[ii],timeArray[jj],satAvg[isnd],satAvg[jsnd],tauMn[itnd],tauMn[jtnd],tauVar[itnd],tauVar[jtnd],rhoLntau[ijtnd]);
				satCovEP[jjii] = satCovEP[iijj];
				aof<<jj<<"\t"<<ii<<"\t"<<satCovEP[iijj]<<"\t"<<rhoLntau[ijtnd]<<endl;
			}
		}
	}
	aof.close();
	return satCovEP;
}

double Saturation::calcSatCovQuad(double time1, double time2, double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double rhoLntau12)
{
	double lntauMn1,lntauMn2,lntauVar1,lntauVar2;
	double sDRange = 3.5, intRange = 0.1;
	int qPoint = 3;
	double* rloc, * weight;
	GaussianQuadrature(qPoint,rloc,weight);

	int lenn = 1;
	transform1(lenn,&tauMn1,&tauVar1,&lntauMn1,&lntauVar1);
	transform1(lenn,&tauMn2,&tauVar2,&lntauMn2,&lntauVar2);
	double lntauMin1, lntauMax1, lntauMin2, lntauMax2, lntauFront;
	double lPoint1,rPoint1,lPoint2,rPoint2;
	lntauMin1 = lPoint1 = lntauMn1 - (sDRange*sqrt(lntauVar1));
	lntauMax1 = lntauMn1 + (sDRange*sqrt(lntauVar1));
	lntauMin2 = lntauMn2 - (sDRange*sqrt(lntauVar2));
	lntauMax2 = lntauMn2 + (sDRange*sqrt(lntauVar2));

	lntauFront = log(fpwStar*t);
	
	double css = 0;
	//if((rhoLntau12<1)&&(rhoLntau12>0.98)){rhoLntau12 = 0.98;}
	while(true)
	{
		rPoint1 = lPoint1+(intRange*sqrt(lntauVar1));
		if(rPoint1>lntauMax1){rPoint1 = lntauMax1;}
		if((lntauFront<rPoint1)&&(lntauFront>lPoint1)){rPoint1 = lntauFront;}
		for(int iq1=0;iq1<qPoint;iq1++)
		{
			double lntauVal1 = ((rPoint1-lPoint1)/2*rloc[iq1])+((rPoint1+lPoint1)/2);
			double tauVal1 = exp(lntauVal1);
			double sVal1;
			t = time1;
			calcSatDet(1,&tauVal1,&sVal1);
			lPoint2 = lntauMin2;
			if(rhoLntau12<1)
			{
				while(true)
				{
					rPoint2 = lPoint2+(intRange*sqrt(lntauVar2));
					if(rPoint2>lntauMax2){rPoint2 = lntauMax2;}
					if((lntauFront<rPoint2)&&(lntauFront>lPoint2)){rPoint2 = lntauFront;}
					for(int iq2=0;iq2<qPoint;iq2++)
					{
						double lntauVal2 = ((rPoint2-lPoint2)/2*rloc[iq2])+((rPoint2+lPoint2)/2);
						double tauVal2 = exp(lntauVal2);
						double sVal2;
						t = time2;
						calcSatDet(1,&tauVal2,&sVal2);
						css += (rPoint1-lPoint1)/2*weight[iq1]*(rPoint2-lPoint2)/2*weight[iq2]*
							(sVal1-satAvg1)*(sVal2-satAvg2)*pdfGauss2(lntauVal1,lntauMn1,lntauVar1,
																		lntauVal2,lntauMn2,lntauVar2,
																		rhoLntau12);
						//css += (sw1-satAvg1)*(sw2-satAvg2)*pdfGauss2(ltau1,tauMnLog1,tauVarLog1,
						//									 ltau2,tauMnLog2,tauVarLog2, 
						//									 rhoLntau12)*dltau1*dltau2;
					}
					if(rPoint2==lntauMax2){break;}
					lPoint2 = rPoint2;
				}
			}
			else
			{
				int sign;
				if(lntauVal1>lntauMn1){sign = 1;}
				else{sign = -1;}
				double lntauVal2 = lntauMn2+(sign*sqrt((lntauVar2/lntauVar1*((lntauVal1-lntauMn1)*(lntauVal1-lntauMn1)))+(2*lntauVar2*log(sqrt(lntauVar1/lntauVar2)))));
				double tauVal2 = exp(lntauVal2);
				double sVal2;
				t = time2;
				calcSatDet(1,&tauVal2,&sVal2);
				css += (rPoint1-lPoint1)/2*weight[iq1]*(sVal1-satAvg1)*(sVal2-satAvg2)*pdfGauss(lntauVal1,lntauMn1,lntauVar1);
			}
		}
		if(rPoint1==lntauMax1){break;}
		lPoint1 = rPoint1;
	}
	delete[] rloc;
	delete[] weight;
	return css;
}

// --- function to calculate fw and fpw at a given sw -----
double Saturation::fw(double sw)
{
  return visR*(sw-Swc)*(sw-Swc)
    /(visR*(sw-Swc)*(sw-Swc)+(1-sw-Sor)*(1-sw-Sor));
}

double Saturation::fpw(double sw)   // dfw_dsw
{
  return (2*visR*deltaS*(sw-Swc)*(1-sw-Sor))
    /pow(visR*(sw-Swc)*(sw-Swc)+(1-sw-Sor)*(1-sw-Sor), 2);
}


// --- calculate the root of (fpw-a=0) ----
// --- a should between 0 and fpwStar, else no right solution -----
double Saturation::root(double a, double eps)
{
  double x0 = 1-Sor, x1 = sStar, xm;

  while(true) {
    xm = (x0+x1)/2.0;
    (fpw(xm)-a > 0)? x1=xm : x0=xm;
    if(fabs(x1-x0)<eps) break;
  }

  return xm;
}

// --- tau's pdf -------------------------
double Saturation::pdfForTau(double ta, double tMn, double tVari)
{
  double ltMn, ltVari;
  ltMn = log(tMn) -0.5*log(1+tVari/tMn/tMn);    // mean of ln(tau)
  ltVari = log(1+tVari/tMn/tMn);                // vari of ln(tau)

  return 1/sqrt(2*PI*ltVari)/ta
    *exp(-0.5/ltVari*pow(log(ta)-ltMn, 2));
}

void Saturation::output(int i) {
  /*char file1[100], file2[100], file3[100] ;
  sprintf( file1, "satDet_i=%d.out", i );
  sprintf( file2, "satAvg_i=%d.out", i );
  sprintf( file3, "satVar_i=%d.out", i );
  int iln, jnd, ij;
  ofstream os1(file1, ios::out);
  for(iln=0; iln<nLn; iln++) {
      int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
      for(jnd = 0; jnd < len_tmp; jnd++) {
          ij = jnd + iln*nNd;
          os1 << x1C[ij] << " " << etaMnC[ij] << " " << satDet[ij] << endl;
      } 
      os1 << endl; 
  }
  os1.close();
  
  ofstream os2(file2, ios::out);
  for(iln=0; iln<nLn; iln++) {
      int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
      for(jnd = 0; jnd < len_tmp; jnd++) {
          ij = jnd+iln*nNd;
          os2 << x1C[ij] << " " << etaMnC[ij] << " " << satAvg[ij] << endl;
      }
      os2 << endl; 
  }
  os2.close();
  
  ofstream os3(file3, ios::out);
  for(iln=0; iln<nLn; iln++) {
      int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
      for(jnd = 0; jnd < len_tmp; jnd++) {
          ij = jnd+iln*nNd;
          os3 << x1C[ij] << " " << etaMnC[ij] << " " << satVar[ij] << endl;
      } 
      os3 << endl; 
  }
  os3.close();*/
	char file1[100], file2[100], file3[100] ;
	sprintf( file1, "satDet_i=%d.out", i );
	sprintf( file2, "satAvg_i=%d.out", i );
	sprintf( file3, "satVar_i=%d.out", i );
	int pointStart, ij, sij;
	ofstream os1(file1, ios::out);
	for(int istrl=0; istrl<tSize; istrl++) {
		int len_tmp = ptclTrackPtr->getParticle(swStrl[istrl])->getLength();
		if(endPointFlag==0){pointStart = 0;}
		else{pointStart = len_tmp-1;}
		for(int jnd = pointStart; jnd < len_tmp; jnd++) {
			ij = jnd + swStrl[istrl]*nNd;
			sij = (istrl*(((1-endPointFlag)*nNd)+(endPointFlag)))+((1-endPointFlag)*jnd);
			os1 << x1C[ij] << " " << etaMnC[ij] << " " << satDet[sij] << endl;
		} 
		os1 << endl; 
	}
	os1.close();
	  
	ofstream os2(file2, ios::out);
	for(int istrl=0; istrl<tSize; istrl++) {
		int len_tmp = ptclTrackPtr->getParticle(swStrl[istrl])->getLength();	
		if(endPointFlag==0){pointStart = 0;}
		else{pointStart = len_tmp-1;}
		for(int jnd = pointStart; jnd < len_tmp; jnd++) {
			ij = jnd + swStrl[istrl]*nNd;
			sij = (istrl*(((1-endPointFlag)*nNd)+(endPointFlag)))+((1-endPointFlag)*jnd);
			os2 << x1C[ij] << " " << etaMnC[ij] << " " << satAvg[sij] << endl;
		}
		os2 << endl; 
	}
	os2.close();
	  
	ofstream os3(file3, ios::out);
	for(int istrl=0; istrl<tSize; istrl++) {
		int len_tmp = ptclTrackPtr->getParticle(swStrl[istrl])->getLength();
		if(endPointFlag==0){pointStart = 0;}
		else{pointStart = len_tmp-1;}
		for(int jnd = pointStart; jnd < len_tmp; jnd++) {
			ij = jnd + swStrl[istrl]*nNd;
			sij = (istrl*(((1-endPointFlag)*nNd)+(endPointFlag)))+((1-endPointFlag)*jnd);
			os3 << x1C[ij] << " " << etaMnC[ij] << " " << satVar[sij] << endl;
		} 
		os3 << endl; 
	}
	os3.close();
}

/*// --- calculate the stochastic solution -----------
void Saturation::calcSatAvgVar() {
  int iln, jnd, ij, it;
  double dtau = fpwStar * t/nS, tmp;  

  // --- calculate saturation profile (tau, s) -----------
  for(it = 0; it < nS; it++) 
      tau[it] = (it + 0.5)*dtau; // always use half point
  calcSatDet(nS, tau, s);
  
  for(it = 0; it < nS; it++) 
      s[it] -= Swc;  // now it store Sw-Swc   

  // --- calculate saturation mean and variance ----------
  for(iln = 0; iln < nLn; iln++) {
      int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
      for(jnd = 0; jnd < len_tmp; jnd++) {
          ij = jnd + iln * nNd;
          if(tauMn[ij]/tauMn[nNd-1+iln*nNd] < 1.0E-2 || 
			         tauVar[ij] < 1.0E-10 ) {
             satAvg[ij] = satDet[ij];   // quite deterministic
             satVar[ij] = 0.0;
          } else {                      // do the integration
             satAvg[ij] = satVar[ij] = 0.0;
             for(it = 0; it < nS; it++) {    
				 tmp = pdfForTau(tau[it], tauMn[ij], tauVar[ij])*s[it]*dtau;
                 satAvg[ij]   += tmp;
                 satVar[ij] += tmp*s[it];
             }
             // --- move back the Swc (because I use s is actually Sw-Swc) ---
             satAvg[ij] += Swc;         
             satVar[ij] += (2*Swc*satAvg[ij]-Swc*Swc)-(satAvg[ij]*satAvg[ij]);
          }
      }
  }
}*/
/*void Saturation::calcSatAvgVar2() {
  int iln, jnd, ij, it;

  // --- calculate saturation mean and variance ----------
  for(iln = 0; iln < nLn; iln++) {
      int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
      for(jnd = 0; jnd < len_tmp; jnd++) {
            ij = jnd + iln * nNd;
			if(tauMn[ij]==0) //tauMn == 0 ==> tauVar == 0
			{
				satAvg[ij] = satDet[ij];   // quite deterministic
				satVar[ij] = 0.0;
			}
			else
			{
				// do the integration
				double lntauMn,lntauVar;
				double lntauMin,lntauMax,dlntau;
				int lenn = 1;
				transform1(lenn,&tauMn[ij],&tauVar[ij],&lntauMn,&lntauVar);
				lntauMin = lntauMn - (3.5*sqrt(lntauVar));
				lntauMax = lntauMn + (3.5*sqrt(lntauVar));
				dlntau = (lntauMax-lntauMin)/nS;
				satAvg[ij] = satVar[ij] = 0.0;
				for(it = 0; it < nS; it++) {
					double lntauVal = (it+0.5)*dlntau+lntauMin;
					double tauVal = exp(lntauVal);
					double sVal;
					calcSatDet(1,&tauVal,&sVal);
					double tmp = pdfGauss(lntauVal,lntauMn,lntauVar)*sVal*dlntau;
					satAvg[ij]   += tmp;
					satVar[ij] += tmp*sVal;
				}    
				satVar[ij] -= (satAvg[ij]*satAvg[ij]);
			}
      }
  }
}*/
/*void Saturation::calcSatAvgVarQuad(){
	int iln, jnd, ij, iq;
	double sDRange = 3.5, intRange = 1;
	int qPoint = 3;
	double* rloc, * weight;
	GaussianQuadrature(qPoint,rloc,weight);

	// --- calculate saturation mean and variance ----------
	for(iln = 0; iln < nLn; iln++) {
		int len_tmp = ptclTrackPtr->getParticle(iln)->getLength();	  
		for(jnd = 0; jnd < len_tmp; jnd++) {
				ij = jnd + iln * nNd;
				if(tauMn[ij]==0) //tauMn == 0 ==> tauVar == 0
				{
					satAvg[ij] = satDet[ij];   // quite deterministic
					satVar[ij] = 0.0;
				}
				else
				{
					// do the integration
					double lntauMn,lntauVar,sVal;
					double lntauMin,lntauMax,lntauFront,lPoint,rPoint;
					int lenn = 1;
					transform1(lenn,&tauMn[ij],&tauVar[ij],&lntauMn,&lntauVar);
					lntauMin = lPoint = lntauMn - (sDRange*sqrt(lntauVar));
					lntauMax = lntauMn + (sDRange*sqrt(lntauVar));
					lntauFront = log(fpwStar*t);
					satAvg[ij] = satVar[ij] = 0.0;
					while(true)
					{
						rPoint = lPoint+(intRange*sqrt(lntauVar));
						if(rPoint>lntauMax){rPoint = lntauMax;}
						if((lntauFront<rPoint)&&(lntauFront>lPoint)){rPoint = lntauFront;}
						for(iq = 0; iq < qPoint; iq++) {
							double lntauVal = ((rPoint-lPoint)/2*rloc[iq])+((rPoint+lPoint)/2);
							double tauVal = exp(lntauVal);
							double sVal;
							calcSatDet(1,&tauVal,&sVal);
							double tmp = (rPoint-lPoint)/2*weight[iq]*sVal*pdfGauss(lntauVal,lntauMn,lntauVar);
							satAvg[ij]   += tmp;
							satVar[ij] += tmp*sVal;
						}
						if(rPoint==lntauMax){break;};
						lPoint = rPoint;
					}
					satVar[ij] -= (satAvg[ij]*satAvg[ij]);
				}
		}
	}
	delete[] rloc;
	delete[] weight;
}
*/
/*double Saturation::calcSatCov(double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double rhoLntau12)
{
	double tauMnLog1,tauMnLog2,tauVarLog1,tauVarLog2;
	
	int lenn = 1;
	transform1(lenn,&tauMn1,&tauVar1,&tauMnLog1,&tauVarLog1);
	transform1(lenn,&tauMn2,&tauVar2,&tauMnLog2,&tauVarLog2);
	double dtau = fpwStar * t/nS;
	
	double css = 0;
	for(int ii=0;ii<nS;ii++)
	{
		double tau1 = (ii+0.5)*dtau;
		double ltau1 = log(tau1);
		double dltau1 = dtau/tau1;
		double sw1;
		calcSatDet(1,&tau1,&sw1);
		for(int jj=0;jj<nS;jj++)
		{
			double tau2 = (jj+0.5)*dtau;
			double ltau2 = log(tau2);
			double dltau2 = dtau/tau2;
			double sw2;
			calcSatDet(1,&tau2,&sw2);
			css += (sw1-Swc)*(sw2-Swc)*pdfGauss2(ltau1,tauMnLog1,tauVarLog1,
												 ltau2,tauMnLog2,tauVarLog2, 
												 rhoLntau12)*dltau1*dltau2;
		}
	}
	css += (satAvg1*Swc) + (satAvg2*Swc) - (Swc*Swc) - (satAvg1*satAvg2);
	return css;
}*/
/*double Saturation::calcSatCov2(double satAvg1, double satAvg2, double tauMn1, double tauMn2, double tauVar1, double tauVar2, double rhoLntau12)
{
	double tauMnLog1,tauMnLog2,tauVarLog1,tauVarLog2;
	
	int lenn = 1;
	transform1(lenn,&tauMn1,&tauVar1,&tauMnLog1,&tauVarLog1);
	transform1(lenn,&tauMn2,&tauVar2,&tauMnLog2,&tauVarLog2);
	double tauLogMin1, tauLogMax1, tauLogMin2, tauLogMax2;
	tauLogMin1 = tauMnLog1 - (3.5*sqrt(tauVarLog1));
	tauLogMax1 = tauMnLog1 + (3.5*sqrt(tauVarLog1));
	tauLogMin2 = tauMnLog2 - (3.5*sqrt(tauVarLog2));
	tauLogMax2 = tauMnLog2 + (3.5*sqrt(tauVarLog2));
	double dltau1 = (tauLogMax1-tauLogMin1)/nS;
	double dltau2 = (tauLogMax2-tauLogMin2)/nS;
	
	double css = 0;
	for(int ii=0;ii<nS;ii++)
	{
		double ltau1 = (ii+0.5)*dltau1+tauLogMin1;
		double tau1 = exp(ltau1);
		double sw1;
		calcSatDet(1,&tau1,&sw1);
		for(int jj=0;jj<nS;jj++)
		{
			double ltau2 = (jj+0.5)*dltau2+tauLogMin2;
			double tau2 = exp(ltau2);
			double sw2;
			calcSatDet(1,&tau2,&sw2);
			css += (sw1-satAvg1)*(sw2-satAvg2)*pdfGauss2(ltau1,tauMnLog1,tauVarLog1,
												 ltau2,tauMnLog2,tauVarLog2, 
												 rhoLntau12)*dltau1*dltau2;
		}
	}
	return css;
}*/

