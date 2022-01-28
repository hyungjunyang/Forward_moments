#include "SwCondiInverse.h"

SwCondiInverse::SwCondiInverse(Domain* domain, Flow* flow)
: domainPtr(domain), flowPtr(flow)
{
	gridPtr = domainPtr->getGrid();
	permPtr = flowPtr->getPermPtr();
	P0e = (P0Eqn*)flowPtr->getP0e();
	CYPe = (CYPEqn*)flowPtr->getCYPe();
	debug = false;
	condiTimePtr = NULL;
	ReadData();
	Initialize();
	//testBackTracking();
	cout<<"SwCondiInverse object created."<<endl;
}

SwCondiInverse::~SwCondiInverse()
{
	if(condiTimePtr != NULL)
	{
		for(int ct=0;ct<totalCondiTime;ct++)
		{
			if(condiTimePtr[ct].numSw>0)
			{
				delete[] condiTimePtr[ct].swPtr;
			}
		}
		delete[] condiTimePtr;
	}
}

void SwCondiInverse::ReadData()
{
	cout<<"Start reading data for SwCondiInverse"<<endl;
	ifstream in_SwCondi;
	in_SwCondi.open("trans_condi.in",ios::in);
	if(in_SwCondi.bad())
	{
		cerr<<"Can not open file, trans_condi.in"<<endl;
		exit(8);
	} 
	else 
	{
		if(debug)
		{
			cout<<"File trans_condi.in was opened successfully!"<<endl;
		}
	}

	Junk(in_SwCondi);
	in_SwCondi >> totalCondiTime;
	condiTimePtr = new CondiTime[totalCondiTime];
	for(int ct=0;ct<totalCondiTime;ct++)
	{
		Junk(in_SwCondi);
		in_SwCondi >> condiTimePtr[ct].ctime >> condiTimePtr[ct].numSw;
		if(condiTimePtr[ct].numSw > 0)
		{
			condiTimePtr[ct].swPtr = new SwData[condiTimePtr[ct].numSw];
		}
		for(int nsw=0;nsw<condiTimePtr[ct].numSw;nsw++)
		{
			Junk(in_SwCondi);
			in_SwCondi >> condiTimePtr[ct].swPtr[nsw].i_ >> condiTimePtr[ct].swPtr[nsw].j_
					   >> condiTimePtr[ct].swPtr[nsw].k_ >> condiTimePtr[ct].swPtr[nsw].sw;
			// --- return back to "0" based index ------
			condiTimePtr[ct].swPtr[nsw].i_--;
			condiTimePtr[ct].swPtr[nsw].j_--;
			condiTimePtr[ct].swPtr[nsw].k_--;
		}
	}
	Junk(in_SwCondi);
	in_SwCondi >> maxIter >> etol;
	cout<<"End reading data for SwCondiInverse"<<endl;
}

void SwCondiInverse::Initialize()
{
	Swc     = domainPtr->getSwc();
	Sor     = domainPtr->getSor();
	viscosR = domainPtr->getViscosR();
	nS = 1000;
	blSolnPtr = new BLSolution(nS, Swc,  Sor, viscosR);

	covSY = NULL;
}

void SwCondiInverse::SequentialConditioning(char *directory, ifstream &is_flow)
{
	cout<<"Start Saturation Conditioning Sequentially"<<endl;
	
	//Sequential Time Loop
	for(int ii=0;ii<totalCondiTime;ii++)
	{
		krigPtr = new KrigWeigh(gridPtr,permPtr,condiTimePtr[ii].numSw);
		//Iteration Loop
		char satFileName[100];
		sprintf(satFileName,"SatMean_Time=%.3f[satinv].out",condiTimePtr[ii].ctime);
		ofstream sof(satFileName, ios::out);
		sof<<"[Time,Iteration#,Particle#] Sw(Measured) Sw(Predicted) Diff Sw_Variance"<<endl<<endl;
		int iter = 0;
		while(true)
		{
			//Get Perm X Y Z
			/*ofstream pxy("PermXYZ.out", ios::out);
			pxy<<"##### Perm X Y #####"<<endl;
			for(int jjj=0;jjj<gridPtr->getPermNy();jjj++)
			{
				for(int iii=0;iii<gridPtr->getPermNx();iii++)
				{
					int iijj = gridPtr->getPermIndex(iii,jjj,0);
					pxy<<iijj<<"   "<<gridPtr->getPermX(iii)<<"   "<<gridPtr->getPermY(jjj)<<endl;
				}
			}
			pxy.close();
			break;*/

			cout<<"Iteration "<<iter<<endl;
			//Lauch paticles from all measurement locations, and calc traveltime moments
			condiTimePtr[ii].numSw = 2;		//For Sequential Observation & Multiple Observation
			LaunchParticles(ii,iter);
			//Solve for saturation moments: mean, var, cov
			CalcSaturationMoments(ii,iter,sof);
			condiTimePtr[ii].numSw = 1;		//For Sequential Observation & Multiple Observation
			//Check whether the conditioning is convert
			if((MaxAbs(condiTimePtr[ii].numSw,satErr)<=etol)||(iter==maxIter))
			{
				DeallocateIter();
				break;
			}
			//Calculate covariance of saturation and log-permeability
			cout<<endl<<endl<<"##### CSY #####"<<endl;
			CalcCovSYVector(ii);
			//Calculate conditioning weight (Lamda)
			CalcLamdaVector(ii);
			//Update permeability moments
			UpdatePermeabilityMoments(ii,iter);
			/*DeallocateIter();
			break;*/
			//Resolve flow problem
			flowPtr->solve(directory,is_flow);
			flowPtr->writeToFile(directory);
			flowPtr->writeToDomain(domainPtr);
			//Dallocate HEAP
			DeallocateIter();
			iter++;
		}
		delete krigPtr;
		sof.close();
	}
}

void SwCondiInverse::LaunchParticles(int tInd, int iter)
{
	pointPtr = new Point[condiTimePtr[tInd].numSw];
	for(int jj=0;jj<condiTimePtr[tInd].numSw;jj++)
	{
		int i_c = condiTimePtr[tInd].swPtr[jj].i_;
		int j_c = condiTimePtr[tInd].swPtr[jj].j_;
		int k_c = condiTimePtr[tInd].swPtr[jj].k_;
		double x_c = gridPtr->getX(i_c);
		double y_c = gridPtr->getY(j_c);
		double z_c = gridPtr->getZ(k_c);
		z_c = 0; //Set to 0 to be consistent
		Point tempPtr(x_c,y_c,z_c,i_c,j_c,k_c);
		pointPtr[jj] = tempPtr;
	}
	bTrackPtr = new ParticleTrack(domainPtr,condiTimePtr[tInd].numSw,pointPtr);
	bTrackPtr->goTrackBackward();
	bTrackPtr->calcTrvlTimeMomentsNoWell();
	bTrackPtr->printTravelTime();
	// Display particle info
	char tauFileName[100];
	sprintf(tauFileName,"TravelTime_Time=%.3f_Iter=%d[satinv].out",condiTimePtr[tInd].ctime,iter);
	ofstream tof(tauFileName, ios::out);
	tof<<"xp yp distance tau_avg0tran tau_avg2tran tau_vartran tau_avg0carte tau_avg2carte tau_varcarte  tau_avg0 tau_avg2 tau_var"<<endl<<endl;
	for(int jj=0;jj<condiTimePtr[tInd].numSw;jj++)
	{
		bTrackPtr->getParticle(jj)->PrintPathline(tof);
	}
	tof.close();
	cout<<"Time "<<condiTimePtr[tInd].ctime<<endl;
	DisplayParticleInfo(condiTimePtr[tInd].numSw,bTrackPtr);
}

void SwCondiInverse::CalcSaturationMoments(int tInd, int iter, ofstream &sof)
{
	//Compute the moments
	satPtr = new Saturation(bTrackPtr, viscosR, Swc, Sor, nS );
	int obstInd;								//For Sequential Observation
	(tInd==0) ? obstInd=1 : obstInd=0;			//For Sequential Observation
	cout<<obstInd<<endl;						//For Sequential Observation
	satPtr->solve(condiTimePtr[obstInd].ctime); //For Sequential Observation
	satPtr->output((tInd+1)*10000+iter);		//For Sequential Observation
	satPtr->solve(condiTimePtr[tInd].ctime);
	satPtr->output((tInd+1)*1000+iter);
	satPtr->calcSatCovEndPoints();
	satCovEP = satPtr->getSatCovEP();

	//Display covariance
	/*cout<<endl<<"Css: ";
	for(int ssind=0;ssind<(condiTimePtr[tInd].numSw*condiTimePtr[tInd].numSw);ssind++)
	{
		cout<<satCovEP[ssind]<<" ";
	}
	cout<<endl<<endl;*/

	//Extract saturation mean
	satAvgEP = new double[condiTimePtr[tInd].numSw];
	satMeas = new double[condiTimePtr[tInd].numSw];
	satErr = new double[condiTimePtr[tInd].numSw];
	for(int ii=0;ii<condiTimePtr[tInd].numSw;ii++)
	{
		satAvgEP[ii] = satPtr->getSatAvg(ii+1,bTrackPtr->getParticle(ii)->getLength());
		satMeas[ii] = condiTimePtr[tInd].swPtr[ii].sw;
		satErr[ii] = satMeas[ii]-satAvgEP[ii];
		cout<<satErr[ii]<<"  ";

		// Output to file
		sof<<"["<<condiTimePtr[tInd].ctime<<","<<iter<<","<<ii<<"] ";
		sof<<satMeas[ii]<<"  "<<satAvgEP[ii]<<"  "<<satErr[ii];
		sof<<"  "<<satCovEP[ii*condiTimePtr[tInd].numSw+ii]<<endl;

	}
	sof<<endl;
	cout<<endl;

	//Display mean
	cout<<endl<<"Sw Mean: ";
	for(int sind=0;sind<condiTimePtr[tInd].numSw;sind++)
	{
		cout<<satAvgEP[sind]<<" ";
	}
	cout<<endl<<endl;
}

void SwCondiInverse::DisplayParticleInfo(int npoint, ParticleTrack* pTrackPtr)
{
	for(int jj=0;jj<npoint;jj++)
		{
			vector<Point> ptraj = pTrackPtr->getParticle(jj)->getVecTrajectory();
			vector<double> ptime = pTrackPtr->getParticle(jj)->getVecTimeFlight();
			int psize = ptraj.size();
			cout<<"POINT "<<jj<<": "<<endl;
			cout<<"Trajectory"<<endl;
			for(int kk=0;kk<psize;kk++)
			{
				Point ppoint = ptraj.back();
				cout<<ppoint.getI()<<" "<<ppoint.getJ()<<" "<<ppoint.getK()<<" "
					<<ppoint.getX()<<" "<<ppoint.getY()<<" "<<ppoint.getZ();
				ptraj.pop_back();
				cout<<endl;
			}
			cout<<endl<<"Time_Flight"<<endl;
			for(int kk=0;kk<psize;kk++)
			{
				cout<<ptime.back()<<" "<<pTrackPtr->getParticle(jj)->getTrvlTimeAvg(psize-1-kk);
				cout<<" "<<pTrackPtr->getParticle(jj)->getTrvlTimeVar(psize-1-kk);
				ptime.pop_back();
				cout<<endl;
			}
		}

}

double SwCondiInverse::MaxAbs(int asize, double* arr)
{
	double maxValue = fabs(arr[0]);
	for(int ii=1;ii<asize;ii++)
	{
		if(fabs(arr[ii])>maxValue){maxValue = fabs(arr[ii]);}
	}
	return maxValue;
}

void SwCondiInverse::DeallocateIter()
{
	delete[] pointPtr;
	delete bTrackPtr;
	delete satPtr;
	delete[] satAvgEP;
	delete[] satMeas;
	delete[] satErr;
	if(covSY!=NULL)
	{
		delete[] covSY;
		covSY = NULL;
		delete[] lamda;
		lamda = NULL;
	}
}

void SwCondiInverse::CalcCovSYVector(int tInd)
{
	covSY = new double[condiTimePtr[tInd].numSw*gridPtr->getNumPermNode()];
	ofstream csy("C_SY.out", ios::out);
	csy<<"##### CSY #####"<<endl;

	for(int sInd=0;sInd<condiTimePtr[tInd].numSw;sInd++)
	{
		CollectDataS(sInd);
		for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
		{
			double st = clock();
			for(int pInd=0;pInd<bTrackPtr->getParticle(sInd)->getLength();pInd++)
			{
				CalcCov_YVxi_YVeta_YEta(sInd,yInd,pInd);
			}
			int syInd = sInd + (yInd*condiTimePtr[tInd].numSw);
			covSY[syInd] = CalcCovSY(sInd,yInd);
			cout<<sInd<<"   "<<yInd<<"   "<<covSY[syInd];
			csy<<sInd<<"   "<<yInd<<"   "<<covSY[syInd];
			double et = clock();
			cout<<"   // Time: "<<et-st<<endl;
			csy<<"   // Time: "<<et-st<<endl;
		}
		DeallocateDataS();
	}
	csy.close();
}

void SwCondiInverse::CollectDataS(int sInd)
{
	int nPoint = bTrackPtr->getParticle(sInd)->getLength();
	ds = new double[nPoint];
	dvdeta = new double[nPoint];
	costh = new double[nPoint];sinth = new double[nPoint];
	x_ = new double [nPoint];y_ = new double [nPoint];z_ = new double [nPoint];
	i_ = new int[nPoint];j_ = new int[nPoint];k_ = new int[nPoint];
	xlb_ = new double [nPoint];ylb_ = new double [nPoint];zlb_ = new double [nPoint];
	xrt_ = new double [nPoint];yrt_ = new double [nPoint];zrt_ = new double [nPoint];
	edge = new int[nPoint];
	covYVxiU = new double[nPoint];covYVxiD = new double[nPoint];
	covYVetaU = new double[nPoint];covYVetaD = new double[nPoint];
	covYEta = new double[nPoint];

	for(int ii=0;ii<nPoint;ii++)
	{
		ds[ii] = bTrackPtr->getParticle(sInd)->getDs(ii);
		dvdeta[ii] = bTrackPtr->getParticle(sInd)->getDvDeta(ii);
		costh[ii] = bTrackPtr->getParticle(sInd)->getCosth(ii);
		sinth[ii] = bTrackPtr->getParticle(sInd)->getSinth(ii);
		FindEdge(sInd,ii);
	}
}

void SwCondiInverse::DeallocateDataS()
{
	delete[] ds;
	delete[] dvdeta;
	delete[] costh;delete[] sinth;
	delete[] x_;delete[] y_;delete[] z_;
	delete[] i_;delete[] j_;delete[] k_;
	delete[] xlb_;delete[] ylb_;delete[] zlb_;
	delete[] xrt_;delete[] yrt_;delete[] zrt_;
	delete[] edge;
	delete[] covYVxiU;delete[] covYVxiD;
	delete[] covYVetaU;delete[] covYVetaD;
	delete[] covYEta;
}

//   5_3_7
//  0| 8 |1
//   4_2_6  edge value (W=0,E=1,S=2,N=3,SW=4,NW=5,SE=6,NE=7)
void SwCondiInverse::FindEdge(int sInd, int pInd)
{
	double xzero = 1e-6;

	x_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getX();
	y_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getY();
	z_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getZ();
	i_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getI();
	j_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getJ();
	k_[pInd] = bTrackPtr->getParticle(sInd)->getTrajectory(pInd).getK();

	xlb_[pInd] = domainPtr->getBlock(i_[pInd],j_[pInd])->getBlockLbPnt().getX();
	ylb_[pInd] = domainPtr->getBlock(i_[pInd],j_[pInd])->getBlockLbPnt().getY();
	xrt_[pInd] = xlb_[pInd] + domainPtr->getBlock(i_[pInd],j_[pInd])->getDx();
	yrt_[pInd] = ylb_[pInd] + domainPtr->getBlock(i_[pInd],j_[pInd])->getDy();

	bool ww,ee,ss,nn; // West, East, South, North
	ww = ee = ss = nn = false;
	if(fabs(x_[pInd]-xlb_[pInd])<xzero){ww = true;}
	if(fabs(x_[pInd]-xrt_[pInd])<xzero){ee = true;}
	if(fabs(y_[pInd]-ylb_[pInd])<xzero){ss = true;}
	if(fabs(y_[pInd]-yrt_[pInd])<xzero){nn = true;}

	if(ww)
	{
		if(ss){edge[pInd] = 4;}
		else if(nn){edge[pInd] = 5;}
		else {edge[pInd] = 0;}
	}
	else if(ee)
	{
		if(ss){edge[pInd] = 6;}
		else if(nn){edge[pInd] = 7;}
		else {edge[pInd] = 1;}
	}
	else if(ss){edge[pInd] = 2;}
	else if(nn){edge[pInd] = 3;}
	else{edge[pInd] = 8;}
}

void SwCondiInverse::CalcCov_YVxi_YVeta_YEta(int sInd, int yInd, int pInd)
{
	double covYVxu,covYVxd,covYVyu,covYVyd; // u = upstream, d = downstream
	double covYVxlu,covYVxld,covYVxru,covYVxrd,covYVybu,covYVybd,covYVytu,covYVytd;
	// l = left, r = right, b = bottom, t = top
	dpdx = P0e->getDP0DXi(0,0);
	dpdy = P0e->getDP0DXi(1,0);

	switch(edge[pInd])
	{
		case 0:
			covYVxlu = covYVxld = covYVxru = covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = NULL;}
			else if(x_[pInd-1]<x_[pInd])
			{
				covYVybu = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
				covYVytu = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
			}
			else
			{
				covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = NULL;}
			else if(x_[pInd+1]>x_[pInd])
			{
				covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			else
			{
				covYVybd = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
				covYVytd = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
			}
			break;
		case 1:
			covYVxlu = covYVxld = covYVxru = covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = NULL;}
			else if(x_[pInd-1]<x_[pInd])
			{
				covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			else
			{
				covYVybu = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				covYVytu = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = NULL;}
			else if(x_[pInd+1]>x_[pInd])
			{
				covYVybd = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				covYVytd = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
			}
			else
			{
				covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			break;
		case 2:
			covYVybu = covYVybd = covYVytu = covYVytd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = NULL;}
			else if(y_[pInd-1]<y_[pInd])
			{
				covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
				covYVxru = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
			}
			else
			{
				covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = NULL;}
			else if(y_[pInd+1]>y_[pInd])
			{
				covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
			}
			else
			{
				covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
				covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
			}
			break;
		case 3:
			covYVybu = covYVybd = covYVytu = covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = NULL;}
			else if(y_[pInd-1]<y_[pInd])
			{
				covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
			}
			else
			{
				covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				covYVxru = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = NULL;}
			else if(y_[pInd+1]>y_[pInd])
			{
				covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
			}
			else
			{
				covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
			}
			break;
		case 4:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 5:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 6:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 7:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else{cerr<<"Particles are align and pass thru a corner, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 8:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
				covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
			}
			break;
		default:
			cerr<<"Point "<<pInd<<": There is error in edge[]"<<endl<<endl;
	}

	if((covYVxlu!=NULL)&&(covYVxru!=NULL))
	{
		covYVxu = ((covYVxlu*(xrt_[pInd]-x_[pInd]))+(covYVxru*(x_[pInd]-xlb_[pInd])))/(xrt_[pInd]-xlb_[pInd]);
	}else{covYVxu = NULL;}
	if((covYVxld!=NULL)&&(covYVxrd!=NULL))
	{
		covYVxd = ((covYVxld*(xrt_[pInd]-x_[pInd]))+(covYVxrd*(x_[pInd]-xlb_[pInd])))/(xrt_[pInd]-xlb_[pInd]);
	}else{covYVxd = NULL;}
	if((covYVybu!=NULL)&&(covYVytu!=NULL))
	{
		covYVyu = ((covYVybu*(yrt_[pInd]-y_[pInd]))+(covYVytu*(y_[pInd]-ylb_[pInd])))/(yrt_[pInd]-ylb_[pInd]);
	}else{covYVyu = NULL;}
	if((covYVybd!=NULL)&&(covYVytd!=NULL))
	{
		covYVyd = ((covYVybd*(yrt_[pInd]-y_[pInd]))+(covYVytd*(y_[pInd]-ylb_[pInd])))/(yrt_[pInd]-ylb_[pInd]);
	}else{covYVyd = NULL;}
	if(pInd == 0){covYVxiU[pInd] = covYVetaU[pInd] = NULL;}
	else
	{
		covYVxiU[pInd] = (covYVxu*costh[pInd])+(covYVyu*sinth[pInd]);
		covYVetaU[pInd] = -(covYVxu*sinth[pInd])+(covYVyu*costh[pInd]);
	}
	if(pInd == bTrackPtr->getParticle(sInd)->getLength()-1)
	{covYVxiD[pInd] = covYVetaD[pInd] = NULL;}
	else
	{
		covYVxiD[pInd] = (covYVxd*costh[pInd])+(covYVyd*sinth[pInd]);
		covYVetaD[pInd] = -(covYVxd*sinth[pInd])+(covYVyd*costh[pInd]);
	}
	if(pInd == 0){covYEta[pInd] = 0;}
	else
	{
		double vxiavg_m1 = bTrackPtr->getParticle(sInd)->getVAvg(pInd-1);
		double vxiavg = bTrackPtr->getParticle(sInd)->getVAvg(pInd);
		covYEta[pInd] = covYEta[pInd-1] + (covYVetaD[pInd-1]/vxiavg_m1*ds[pInd]/2) + (covYVetaU[pInd]/vxiavg*ds[pInd]/2);
	}
}

double SwCondiInverse::CalcCovYVx(int i1, int j1, int k1, int i2, int j2, int k2, int yInd)
{
	double perm_interface[3];
	int permInd, ijkpermtemp_;
	double dpdxtemp, cyytemp, dcpydxtemp;
	int ijk1, ijk2;

	dpdxtemp = dcpydxtemp = 0; // Will not be changed if the interface is at a border (no flow)
	ijk1 = gridPtr->getIndex(i1,j1,k1);
	ijk2 = gridPtr->getIndex(i2,j2,k2);

	if(i1>=0)
	{
		permPtr->Perm_x(i1,j1,k1,perm_interface);
		permInd = 2;
		if(i1<(gridPtr->getNx()-1))
		{
			dpdxtemp = dpdx[gridPtr->getIndex(i1,j1,k1)];
			dcpydxtemp = (CYPe->getCYP(ijk2,yInd)-CYPe->getCYP(ijk1,yInd))/gridPtr->getDx(i1);
		}
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i1)+1,gridPtr->toJPerm(j1),gridPtr->toKPerm(k1));
	}
	else
	{
		permPtr->Perm_x(i2,j2,k2,perm_interface);
		permInd = 0;
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i2)-1,gridPtr->toJPerm(j2),gridPtr->toKPerm(k2));
	}

	cyytemp = permPtr->getCYY(ijkpermtemp_,yInd);
	return(perm_interface[permInd]*((dpdxtemp*cyytemp)+dcpydxtemp));
}

double SwCondiInverse::CalcCovYVy(int i1, int j1, int k1, int i2, int j2, int k2, int yInd)
{
	double perm_interface[3];
	int permInd, ijkpermtemp_;
	double dpdytemp, cyytemp, dcpydytemp;
	int ijk1, ijk2;

	dpdytemp = dcpydytemp = 0; // Will not be changed if the interface is at a border (no flow)
	ijk1 = gridPtr->getIndex(i1,j1,k1);
	ijk2 = gridPtr->getIndex(i2,j2,k2);

	if(j1>=0)
	{
		permPtr->Perm_y(i1,j1,k1,perm_interface);
		permInd = 2;
		if(j1<(gridPtr->getNy()-1))
		{
			dpdytemp = dpdy[gridPtr->getIndex(i1,j1,k1)];
			dcpydytemp = (CYPe->getCYP(ijk2,yInd)-CYPe->getCYP(ijk1,yInd))/gridPtr->getDy(i1);
		}
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i1),gridPtr->toJPerm(j1)+1,gridPtr->toKPerm(k1));
	}
	else
	{
		permPtr->Perm_y(i2,j2,k2,perm_interface);
		permInd = 0;
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i2),gridPtr->toJPerm(j2)-1,gridPtr->toKPerm(k2));
	}

	cyytemp = permPtr->getCYY(ijkpermtemp_,yInd);
	return(perm_interface[permInd]*((dpdytemp*cyytemp)+dcpydytemp));
}

double SwCondiInverse::CalcCovSY(int sInd, int yInd)
{
	double savg = satAvgEP[sInd];
	double yavg = permPtr->getYAvg(yInd);
	double yvar = permPtr->getYVar(yInd);
	double tauavg = bTrackPtr->getParticle(sInd)->getLastTravelTime();
	double tauvar = bTrackPtr->getParticle(sInd)->getLastTravelTVar();
	double lntauavg, lntauvar;
	int lenn = 1;
	transform1(lenn,&tauavg,&tauvar,&lntauavg,&lntauvar);
	double covYTau = CalcCovYTau(sInd);
	double covYLnTau = CalcCovYLnTau(yavg,yvar,tauavg,tauvar,covYTau,1000,1e-6);
	double rhoYLnTau = covYLnTau/sqrt(yvar*lntauvar);

	//Double integration
	double ymin = yavg - (3.5*sqrt(yvar));
	double ymax = yavg + (3.5*sqrt(yvar));
	double dy = (ymax-ymin)/nS;
	double lntaumin = lntauavg - (3.5*sqrt(lntauvar));
	double lntaumax = lntauavg + (3.5*sqrt(lntauvar));
	double dlntau = (lntaumax-lntaumin)/nS;

	double covSY = 0;
	for(int yy=0;yy<nS;yy++)
	{
		double yval = (yy+0.5)*dy + ymin;
		for(int tt=0;tt<nS;tt++)
		{
			double lntauval = (tt+0.5)*dlntau + lntaumin;
			double tauval = exp(lntauval);
			double sval;
			satPtr->calcSatDet(1,&tauval,&sval);
			double lntsd = (lntauval-lntauavg)/sqrt(lntauvar);
			double ysd = (yval-yavg)/sqrt(yvar);
			double param = (sval-savg)*(yval-yavg);
			covSY += param*(1/(2*3.1415926535897932*sqrt(lntauvar*yvar*(1-(rhoYLnTau*rhoYLnTau))))*exp(-0.5/(1-(rhoYLnTau*rhoYLnTau))*((lntsd*lntsd)-(2*rhoYLnTau*lntsd*ysd)+(ysd*ysd))))*dy*dlntau;
		}
	}
	return covSY;
}

double SwCondiInverse::CalcCovYTau(int sInd)
{
	double covYTau = 0;
	for(int pInd=1;pInd<bTrackPtr->getParticle(sInd)->getLength();pInd++)
	{
		double vxiavg = bTrackPtr->getParticle(sInd)->getVAvg(pInd-1);
		covYTau += 1/(vxiavg*vxiavg)*(covYVxiD[pInd-1]+(covYEta[pInd-1]*dvdeta[pInd-1]))*ds[pInd]/2;
		vxiavg = bTrackPtr->getParticle(sInd)->getVAvg(pInd);
		covYTau += 1/(vxiavg*vxiavg)*(covYVxiU[pInd]+(covYEta[pInd]*dvdeta[pInd]))*ds[pInd]/2;
	}
	return covYTau;
}

/*double SwCondiInverse::CalcCovYLnTau(double yMean, double yVar, double tauMean, double tauVar, double covYTau, int ns, double eps)
{
	double covYTauMin,covYTauMax,covYLnTauMin,covYLnTauMax,factor;
	double covYTaui = CalcCovYTauByIntegration(yMean,yVar,tauMean,tauVar,covYTau,ns);//Use covYTau as intial guess of covYLnTau
	if(covYTaui<covYTau)
	{
		covYTauMin = covYTaui;
		covYLnTauMin = covYTau;
		if(covYTau>=0){factor = 10;}
		else{factor = 0.1;}
	}
	else if(covYTaui>covYTau)
	{
		covYTauMax = covYTaui;
		covYLnTauMax = covYTau;
		if(covYTau>=0){factor = 0.1;}
		else{factor = 10;}
	}else{return covYTau;}
	double covYTauiter,covYLnTauiter = covYTau;
	while(true)
	{
		covYLnTauiter *= factor;
		covYTauiter = CalcCovYTauByIntegration(yMean,yVar,tauMean,tauVar,covYLnTauiter,ns);
		if((covYTauiter-covYTau)*(covYTaui-covYTau)<0){break;}
	}
	if(covYTauiter<covYTau)
	{
		covYTauMin = covYTauiter;
		covYLnTauMin = covYLnTauiter;
	}
	else if(covYTauiter>covYTau)
	{
		covYTauMax = covYTauiter;
		covYLnTauMax = covYLnTauiter;
	}else{return covYLnTauiter;}
	return(RecCalcCovYLnTau(yMean,yVar,tauMean,tauVar,covYTau,ns,eps,covYTauMin,covYLnTauMin,covYTauMax,covYLnTauMax));
}*/

double SwCondiInverse::CalcCovYLnTau(double yMean, double yVar, double tauMean, double tauVar, double covYTau, int ns, double eps)
{
	//This implementation is to test analytical convertion, replacing the above computation
	double lntauMean = (2*log(tauMean))-(log(pow(tauMean,2)+tauVar)/2);
	double lntauVar = log(1+(tauVar/pow(tauMean,2)));
	return((1/exp((lntauVar/2)+lntauMean)*(covYTau+(yMean*tauMean)))-yMean);
}

/*double SwCondiInverse::RecCalcCovYLnTau(double yMean, double yVar, double tauMean, double tauVar, double covYTau, int ns, double eps, double covYTauMin, double covYLnTauMin, double covYTauMax, double covYLnTauMax)
{
	if(fabs((covYLnTauMax-covYLnTauMin)/covYLnTauMin)<eps){return((covYLnTauMax+covYLnTauMin)/2);}
	//if(fabs((covYTauMax-covYTau)/covYTau)<eps){return(covYLnTauMax);}
	//if(fabs((covYTauMin-covYTau)/covYTau)<eps){return(covYLnTauMin);}
	double covYTauMid,covYLnTauMid;
	covYLnTauMid = (covYLnTauMax-covYLnTauMin)/(covYTauMax-covYTauMin)*(covYTau-covYTauMin)+covYLnTauMin;
	if((covYLnTauMid==covYLnTauMin)||(covYLnTauMid==covYLnTauMax))
	{covYLnTauMid = (covYLnTauMax+covYLnTauMin)/2;}
	covYTauMid = CalcCovYTauByIntegration(yMean,yVar,tauMean,tauVar,covYLnTauMid,ns);
	if(covYTauMid==covYTau){return covYLnTauMid;}
	if((covYTauMid-covYTau)*(covYTauMin-covYTau)<0)
	{
		covYTauMax = covYTauMid;
		covYLnTauMax = covYLnTauMid;
	}
	else
	{
		covYTauMin = covYTauMid;
		covYLnTauMin = covYLnTauMid;
	}
	return(RecCalcCovYLnTau(yMean,yVar,tauMean,tauVar,covYTau,ns,eps,covYTauMin,covYLnTauMin,covYTauMax,covYLnTauMax));
}*/

/*double SwCondiInverse::CalcCovYTauByIntegration(double yMean, double yVar, double tauMean, double tauVar, double covYLnTau, int ns)
{
	double lnTauMean, lnTauVar;
	int lenn = 1;
	transform1(lenn,&tauMean,&tauVar,&lnTauMean,&lnTauVar);
	double rhoYLnTau = covYLnTau/sqrt(lnTauVar * yVar);
	double yMin = yMean - (3.5*sqrt(yVar));
	double yMax = yMean + (3.5*sqrt(yVar));
	double dy = (yMax-yMin)/ns;
	double lnTauMin = lnTauMean-(3.5*sqrt(lnTauVar));
	double lnTauMax = lnTauMean+(3.5*sqrt(lnTauVar));
	double dlntau = (lnTauMax-lnTauMin)/ns;

	double value = 0;
	for(int yInd=0;yInd<ns;yInd++)
	{
		double yVal = (yInd+0.5)*dy + yMin;
		for(int lntInd=0;lntInd<ns;lntInd++)
		{
			double lntauVal = (lntInd+0.5)*dlntau + lnTauMin;
			double xx = (lntauVal-lnTauMean)/sqrt(lnTauVar);
			double yy = (yVal-yMean)/sqrt(yVar);
			double param = (exp(lntauVal)-tauMean)*(yVal-yMean);
			//double param = (exp(lntauVal))*(yVal);
			value += param*(1/(2*3.1415926535897932*sqrt(lnTauVar*yVar*(1-(rhoYLnTau*rhoYLnTau))))*exp(-0.5/(1-(rhoYLnTau*rhoYLnTau))*((xx*xx)-(2*rhoYLnTau*xx*yy)+(yy*yy))))*dy*dlntau;
		}
	}
	//value -= tauMean*yMean;

	return value;
}*/

void SwCondiInverse::CalcLamdaVector(int tInd)
{
	lamda = new double[condiTimePtr[tInd].numSw*gridPtr->getNumPermNode()];
	double* satCovEPInv = new double[condiTimePtr[tInd].numSw*condiTimePtr[tInd].numSw];
	krigPtr->calCovInverse(condiTimePtr[tInd].numSw,satCovEP,satCovEPInv);
	int numSw = condiTimePtr[tInd].numSw;

	for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
	{
		for(int sInd1=0;sInd1<numSw;sInd1++)
		{
			lamda[sInd1+(yInd*numSw)] = 0;
			for(int sInd2=0;sInd2<numSw;sInd2++)
			{
				lamda[sInd1+(yInd*numSw)] += satCovEPInv[sInd2+(sInd1*numSw)]*covSY[sInd2+(yInd*numSw)];
			}
		}
	}
}

void SwCondiInverse::UpdatePermeabilityMoments(int tInd, int iter)
{
	int numY = gridPtr->getNumPermNode();
	int numSw = condiTimePtr[tInd].numSw;
	char ymFileName[1000], yvFileName[1000];
	sprintf(ymFileName,"YMean_Time=%.3f_Iter=%d.out",condiTimePtr[tInd].ctime,iter+1);
	ofstream ym(ymFileName, ios::out);
	ym<<"##### Y Mean #####"<<endl;
	sprintf(yvFileName,"YVar_Time=%.3f_Iter=%d.out",condiTimePtr[tInd].ctime,iter+1);
	ofstream yv(yvFileName, ios::out);
	yv<<"##### Y Variance #####"<<endl;

	for(int yInd1=0;yInd1<numY;yInd1++)
	{
		//Mean update
		double meann = permPtr->getYAvg(yInd1);
		for(int sInd=0;sInd<numSw;sInd++)
		{
			meann += lamda[sInd+(yInd1*numSw)]*satErr[sInd];
		}
		permPtr->setYAvg(yInd1,meann);
		ym<<yInd1<<"   "<<meann<<endl;


		//CYY update
		for(int yInd2=yInd1;yInd2<numY;yInd2++)
		{
			double lamdaCSY = 0;
			for(int sInd=0;sInd<numSw;sInd++)
			{
				lamdaCSY += lamda[sInd+(yInd1*numSw)]*covSY[sInd+(yInd2*numSw)];
			}
			permPtr->setCYY_Cond(yInd1,yInd2,lamdaCSY);
			if(yInd1!=yInd2){permPtr->setCYY_Cond(yInd2,yInd1,lamdaCSY);}
			if(yInd1==yInd2){yv<<yInd1<<"   "<<permPtr->getCYY(yInd1,yInd2)-lamdaCSY<<endl;}
		}
	}
	ym.close();
	yv.close();
}

/*
void SwCondiInverse::testBackTracking()//This is needed to test backtracking and reorder data only
{
	cout<<"Start testing backtracking"<<endl;
	for(int ii=0;ii<totalCondiTime;ii++)
	{
		Point* pointPtr = new Point[condiTimePtr[ii].numSw];
		for(int jj=0;jj<condiTimePtr[ii].numSw;jj++)
		{
			int i_c = condiTimePtr[ii].swPtr[jj].i_;
			int j_c = condiTimePtr[ii].swPtr[jj].j_;
			int k_c = condiTimePtr[ii].swPtr[jj].k_;
			double x_c = gridPtr->getX(i_c);
			double y_c = gridPtr->getY(j_c);
			double z_c = gridPtr->getZ(k_c);
			Point tempPtr(x_c,y_c,z_c,i_c,j_c,k_c);
			pointPtr[jj] = tempPtr;
		}
		ParticleTrack* bTrackPtr = new ParticleTrack(domainPtr,condiTimePtr[ii].numSw,pointPtr);
		bTrackPtr->goTrackBackward();
		Point* fpointPtr = new Point[condiTimePtr[ii].numSw];
		for(int jj=0;jj<condiTimePtr[ii].numSw;jj++)
		{
			fpointPtr[jj] = bTrackPtr->getParticle(jj)->getParticle_Point();
		}
		ParticleTrack* fTrackPtr = new ParticleTrack(domainPtr,condiTimePtr[ii].numSw,fpointPtr);
		fTrackPtr->goTrack();

		for(int jj=0;jj<condiTimePtr[ii].numSw;jj++)
		{
			vector<Point> btraj = bTrackPtr->getParticle(jj)->getVecTrajectory();
			vector<double> btime = bTrackPtr->getParticle(jj)->getVecTimeFlight();
			vector<Point> ftraj = fTrackPtr->getParticle(jj)->getVecTrajectory();
			vector<double> ftime = fTrackPtr->getParticle(jj)->getVecTimeFlight();
			int fsize = ftraj.size();
			int bsize = btraj.size();
			cout<<"Time "<<condiTimePtr[ii].ctime<<" Point "<<jj<<endl;
			cout<<"Trajectory"<<endl;
			//for (int kk=0;kk<(fsize-bsize)-1;kk++)
			//{
			//	ftraj.pop_back();ftime.pop_back();
			//}
			int loopsize = fsize;
			for(int kk=0;kk<loopsize;kk++)
			{
				Point fpoint = ftraj.back();Point bpoint = btraj.back();
				ftraj.pop_back();
				cout<<fpoint.getI()<<" "<<fpoint.getJ()<<" "<<fpoint.getK()<<" "
					<<fpoint.getX()<<" "<<fpoint.getY()<<" "<<fpoint.getZ()<<"     ";
				if(kk<bsize){
					cout<<bpoint.getI()<<" "<<bpoint.getJ()<<" "<<bpoint.getK()<<" "
						<<bpoint.getX()<<" "<<bpoint.getY()<<" "<<bpoint.getZ();
					btraj.pop_back();
				}
				cout<<endl;
			}
			cout<<endl<<"Time_Flight"<<endl;
			for(int kk=0;kk<loopsize;kk++)
			{
				cout<<ftime.back()<<"  ";
				if(kk<bsize){
					cout<<btime.back();
					btime.pop_back();
				}
				cout<<endl;
				ftime.pop_back();
			}
		}

		delete[] pointPtr;
		delete[] fpointPtr;
		delete bTrackPtr;
		delete fTrackPtr;
	}
	cout<<"Finish testing backtracking, continue?(y/y): ";
	char cans;
	cin>>cans;
}
*/