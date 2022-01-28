#include "SwCondiInverse.h"

SwCondiInverse::SwCondiInverse(Domain* domain, Flow* flow)
: domainPtr(domain), flowPtr(flow)
{
	debug = true;
	nS = 1000;
	readDataFlag = 0;

	gridPtr = domainPtr->getGrid();
	permPtr = flowPtr->getPermPtr();
	P0e = (P0Eqn*)flowPtr->getP0e();
	CYPe = (CYPEqn*)flowPtr->getCYPe();
	Swc     = domainPtr->getSwc();
	Sor     = domainPtr->getSor();
	viscosR = domainPtr->getViscosR();
	blSolnPtr = new BLSolution(nS, Swc,  Sor, viscosR);
	
	coordSet = NULL;
	condSet = NULL;
	obvSet = NULL;
	covSY = NULL;
	condSatPtr = NULL;
	obvSatPtr = NULL;
	obvSwStrl = NULL;
	satCovEP = NULL;
	bTrackPtr = NULL;
	
	//testBackTracking();
	//cout<<"SwCondiInverse object created."<<endl;
}

SwCondiInverse::~SwCondiInverse()
{
	delete blSolnPtr;
	if(readDataFlag == 1)
	{
		if(coordSet != NULL)
		{
			Coord* delPtr;
			Coord* nexPtr;
			nexPtr = coordSet;
			for(int idel=0;idel<nCoord;idel++)
			{
				delPtr = nexPtr;
				nexPtr = delPtr->coordPtr;
				delete delPtr;
			}
		}
		if(condSet != NULL)
		{
			delete[] condSet;
			delete[] condTime;
			delete[] swMeas;
			delete[] swMeasErr;
		}
		if(obvSet != NULL)
		{
			delete[] obvSet;
			delete[] obvTime;
		}
	}
}

void SwCondiInverse::SaturationConditioning(char* directory, ifstream &is_flow)
{
	permPtr->wrtPermMomnt("UnCdPerm");
	clockStart = clock();
	ReadData();
	if(condScheme==0||condScheme==2){SequentialConditioning(directory,is_flow);}
	else if(condScheme==1){GlobalConditioning(directory,is_flow);}
	else{cout<<"This choice of conditioning scheme is not implemented as yet!!!!";}
	clockEnd = clock();
	cout<<endl<<"SATURATION CONDITIONING TOTAL RUNTIME: "<<(clockEnd-clockStart)/CLOCKS_PER_SEC<<"  SECONDS"<<endl;
	permPtr->wrtPermMomnt("CondPerm");
}

void SwCondiInverse::ReadData()
{
	cout<<"Start reading data for SwCondiInverse"<<endl;
	ifstream in_SwCondi;
	string iline;
	stringstream istream;
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
	getline(in_SwCondi,iline);
	istream.str(iline);
	istream.clear();
	istream >> condScheme;
	if (!(istream >> dFactor)){
		dFactor = 1;
	}
	if (!(istream >> meas_err_opt)){
		meas_err_opt = 0;
	}

	Junk(in_SwCondi);
	in_SwCondi >> nCond >> nObv;

	Coord coordTemp;
	condSet = new Coord*[nCond];
	condTime = new double[nCond];
	swMeas = new double[nCond];
	swMeasErr = new double[nCond * nCond];
	nCoord = 0;
	for(int icond=0;icond<nCond;icond++)
	{
		Junk(in_SwCondi);
		getline(in_SwCondi,iline);
		istream.str(iline);
		istream.clear();
		istream >> condTime[icond];
		istream >> coordTemp.i_;
		istream >> coordTemp.j_;
		istream >> coordTemp.k_;
		istream >> swMeas[icond];
		if (meas_err_opt > 0){
			for (int ierr=0;ierr<=icond;ierr++){
				istream >> swMeasErr[icond*nCond+ierr];
				swMeasErr[ierr*nCond+icond] = swMeasErr[icond*nCond+ierr];
			}
		}	
		// --- return back to "0" based index ------
		coordTemp.i_--;
		coordTemp.j_--;
		coordTemp.k_--;

		coordTemp.x_ = gridPtr->getX(coordTemp.i_);
		coordTemp.y_ = gridPtr->getY(coordTemp.j_);
		coordTemp.z_ = gridPtr->getZ(coordTemp.k_);
		condSet[icond] = FindCoordinate(coordSet,coordTemp);
	}
	obvSet = new Coord*[nObv];
	obvTime = new double[nObv];
	for(int iobv=0;iobv<nObv;iobv++)
	{
		Junk(in_SwCondi);
		in_SwCondi >> obvTime[iobv] >> coordTemp.i_ >> coordTemp.j_ >> coordTemp.k_;
		// --- return back to "0" based index ------
		coordTemp.i_--;
		coordTemp.j_--;
		coordTemp.k_--;

		coordTemp.x_ = gridPtr->getX(coordTemp.i_);
		coordTemp.y_ = gridPtr->getY(coordTemp.j_);
		coordTemp.z_ = gridPtr->getZ(coordTemp.k_);
		obvSet[iobv] = FindCoordinate(coordSet,coordTemp);
	}
	Junk(in_SwCondi);
	vtol = 9999999; //Set vtol to infinity
	in_SwCondi >> maxIter >> etol >> vtol; //vtol will remain infinity if there is no input
	cout<<"End reading data for SwCondiInverse"<<endl;
	readDataFlag = 1;
	if (debug){
		cout << endl;
		cout << "Conditioning scheme: " << condScheme << endl;
		cout << "Damping factor: " << dFactor << endl;
		cout << "Meas error option: " << meas_err_opt << endl;
		cout << "nCond nObv: " << nCond << " " << nObv << endl;
		if (meas_err_opt > 0){
			cout << "Measurement error:" << endl;
			for (int ierr1=0;ierr1<nCond;ierr1++){
				for (int ierr2=0;ierr2<nCond;ierr2++){
					cout << swMeasErr[ierr1*nCond+ierr2] << " ";
				}
				cout << endl;
			}
		}
	}
}

SwCondiInverse::Coord* SwCondiInverse::FindCoordinate(Coord* &coordPtr, Coord &coordVal)
{
	if(coordPtr==NULL)
	{
		AddCoordinate(coordPtr,coordVal);
		return coordPtr;
	}
	int comp = CompareCoordinate(coordVal,*coordPtr);
	if(comp==-1)
	{
		AddCoordinate(coordPtr,coordVal);
		return coordPtr;
	}
	if(comp==0){return coordPtr;}
	if(comp==1){return FindCoordinate(coordPtr->coordPtr,coordVal);}
}

int SwCondiInverse::CompareCoordinate(Coord coord1, Coord coord2)
{
	if(coord1.i_<coord2.i_){return -1;}
	if(coord1.i_>coord2.i_){return +1;}
	if(coord1.j_<coord2.j_){return -1;}
	if(coord1.j_>coord2.j_){return +1;}
	if(coord1.k_<coord2.k_){return -1;}
	if(coord1.k_>coord2.k_){return +1;}
	if(coord1.x_<coord2.x_){return -1;}
	if(coord1.x_>coord2.x_){return +1;}
	if(coord1.y_<coord2.y_){return -1;}
	if(coord1.y_>coord2.y_){return +1;}
	if(coord1.z_<coord2.z_){return -1;}
	if(coord1.z_>coord2.z_){return +1;}
	return 0;
}

void SwCondiInverse::AddCoordinate(Coord* &coordPtr, Coord coordVal)
{
	Coord* tempPtr = coordPtr;
	coordPtr = new Coord;
	coordPtr->i_ = coordVal.i_;
	coordPtr->j_ = coordVal.j_;
	coordPtr->k_ = coordVal.k_;
	coordPtr->x_ = coordVal.x_;
	coordPtr->y_ = coordVal.y_;
	coordPtr->z_ = coordVal.z_;
	coordPtr->coordPtr = tempPtr;
	nCoord++;
}

void SwCondiInverse::SequentialConditioning(char *directory, ifstream &is_flow)
{
	cout<<"Start Saturation Conditioning Sequentially"<<endl;
	
	//Sequential Time Loop
	int seqLoc = 0;
	int ibatch = 0;
	//for(int ii=0;ii<totalCondiTime;ii++)
	while(true)
	{
		if(seqLoc==nCond){break;}
		int numSw;
		if(condScheme==0){numSw = CountMeasurements(seqLoc);}
		else{numSw = 1;}
		krigPtr = new KrigWeigh(gridPtr,permPtr,numSw);

		// Output File
		char satFileName[100];
		sprintf(satFileName,"CondSatMoments_Time=%.3f[satinv].out",condTime[seqLoc]);
		if(condScheme = 2){sprintf(satFileName,"CondSatMoments_Time=%.3f_%d[satinv].out",condTime[seqLoc],ibatch);}
		ofstream cof(satFileName, ios::out);
		cof<<"[Time,Iteration#,Particle#] Sw(Measured) Sw(Predicted) Diff Sw_Variance"<<endl<<endl;
		sprintf(satFileName,"ObvSatMoments_Time=%.3f[satinv].out",condTime[seqLoc]);
		if(condScheme = 2){sprintf(satFileName,"ObvSatMoments_Time=%.3f_%d[satinv].out",condTime[seqLoc],ibatch);}
		ofstream oof(satFileName, ios::out);
		oof<<"[Time,Iteration#,Particle#] Sw(Predicted) Sw_Variance"<<endl<<endl;

		//Iteration Loop
		int iter = 0;
		while(true)
		{
			cout<<"Iteration "<<iter<<endl;
			//Lauch paticles from all measurement locations, and calc traveltime moments
			//condiTimePtr[ii].numSw = 2;		//For Sequential Observation & Multiple Observation
			condSwStrl = new int[numSw];
			if(nObv>0){obvSwStrl = new int[nObv];}
			if(true) //Observe at every iterations
			//if((nObv>0)&&((seqLoc==0)&&(iter==0))) //Observe only at the first iterations
			  {LaunchParticles(&condSet[seqLoc],numSw,obvSet,nObv,iter);}
			else {LaunchParticles(&condSet[seqLoc],numSw,NULL,0,iter);}
			//Solve for saturation moments: mean, var, cov
			CalcSaturationMoments(&condSet[seqLoc],seqLoc,numSw,cof,iter,1);//(seqLoc,numSw,cof,0,0,oof,iter);//(seqLoc,numSw,iter,sof);
			if(true) //Observe at every iterations
			//if((nObv>0)&&((seqLoc==0)&&(iter==0))) //Observe only at the first iterations
			  {CalcSaturationMoments(obvSet,0,nObv,oof,iter,0);}
			//Check whether the conditioning is convert
			satCovEP = condSatPtr->calcSatCovEndPoints();
			double* satCovErr = new double[numSw];
	                for (int ii=0;ii<numSw;ii++){satCovErr[ii] = 0;}
	                for (int ir=0;ir<numSw;ir++){
	                        if (meas_err_opt == 0){
	                                satCovErr[ir] = satVarEP[ir];
	                        }
	                        else{
	                                for (int ic=0;ic<numSw;ic++){
        	                                satCovErr[ir] = satCovErr[ir]+fabs(satCovEP[ir*numSw+ic]-swMeasErr[(ir+seqLoc)*nCond+(ic+seqLoc)]);
	                                }
	                        }
	                }

			if(((MaxAbs(numSw,satErr)<=etol)&&(MaxAbs(numSw,satCovErr)<=vtol))||(iter==maxIter))
			{
			  delete[] satCovErr;
			  if(nObv>0){
			    delete[] pointPtr;
			    delete bTrackPtr;
			    LaunchParticles(NULL,0,obvSet,nObv,iter);
			    CalcSaturationMoments(obvSet,0,nObv,oof,iter,0);
			  }
			  DeallocateIter();
			  break;
			}
			delete[] satCovErr;
			if (meas_err_opt == 1){
				for (int ierr1=0;ierr1<numSw;ierr1++){
					for (int ierr2=0;ierr2<numSw;ierr2++){
						satCovEP[ierr1*numSw+ierr2] = satCovEP[ierr1*numSw+ierr2]+swMeasErr[(ierr1+seqLoc)*nCond+(ierr2+seqLoc)];
					}
				}
			}
			//Calculate covariance of saturation and log-permeability
			CalcCovSYVector((seqLoc*maxIter)+iter,numSw,seqLoc);
			//Calculate conditioning weight (Lambda)
			CalcLambdaVector(numSw);
			//Update permeability moments
			UpdatePermeabilityMoments(seqLoc,numSw,iter);
			/*DeallocateIter();
			break;*/
			//Resolve flow problem
			flowPtr->deleteVeloPtr(); // To prevent memory leak by calling flowPtr->solve repeatedly
			flowPtr->solve(directory,is_flow);
			flowPtr->writeToFile(directory);
			domainPtr->deleteVeloCovariances(); // To prevent memory leak by calling flowPtr->writeToDomain repeatedly
			flowPtr->writeToDomain(domainPtr);
			//Dallocate HEAP
			DeallocateIter();
			iter++;
		}
		clockEnd = clock();
		delete krigPtr;
		cof<<endl<<"SATURATION CONDITIONING RUNTIME: "<<(clockEnd-clockStart)/CLOCKS_PER_SEC<<"  SECONDS";
		cof.close();
		oof.clear();
		char condYFileName[100];
		sprintf(condYFileName, "CondPerm%d", ibatch);
		permPtr->wrtPermMomnt(condYFileName);
		seqLoc += numSw;
		ibatch++;
	}
}

void SwCondiInverse::GlobalConditioning(char *directory, ifstream &is_flow)
{
	cout<<"Start Saturation Conditioning Globally"<<endl;

	krigPtr = new KrigWeigh(gridPtr,permPtr,nCond);
	// Output files
	char satFileName[100];
	sprintf(satFileName,"CondSatMoments_GlobalScheme[satinv].out");
	ofstream cof(satFileName, ios::out);
	cof<<"[Time,Iteration#,Particle#] Sw(Measured) Sw(Predicted) Diff Sw_Variance"<<endl<<endl;
	sprintf(satFileName,"ObvSatMoments_GlobalScheme[satinv].out");
	ofstream oof(satFileName, ios::out);
	oof<<"[Time,Iteration#,Particle#] Sw(Predicted) Sw_Variance"<<endl<<endl;

	//Iteration
	int iter = 0;
	while(true)
	{
		cout<<"Iteration "<<iter<<endl;
		condSwStrl = new int[nCond];
		if(nObv>0){obvSwStrl = new int[nObv];}
		//Launch particles from all locations at all times
		if((nObv>0)&&(iter==0))
		  {LaunchParticles(condSet,nCond,obvSet,nObv,iter);}
		else {LaunchParticles(condSet,nCond,NULL,0,iter);}
		//Solve for saturation moments: mean, var, cov
		CalcSaturationMoments(condSet,0,nCond,cof,iter,1);
		if((nObv>0)&&(iter==0))
		  {CalcSaturationMoments(obvSet,0,nObv,oof,iter,0);}
		//Check whether the conditioning is convert
		satCovEP = condSatPtr->calcSatCovEndPoints();
		double* satCovErr = new double[nCond];
		for (int ii=0;ii<nCond;ii++){satCovErr[ii] = 0;}
		for (int ir=0;ir<nCond;ir++){
			if (meas_err_opt == 0){
				satCovErr[ir] = satVarEP[ir];
			}
			else{
				for (int ic=0;ic<nCond;ic++){
					satCovErr[ir] = satCovErr[ir]+fabs(satCovEP[ir*nCond+ic]-swMeasErr[ir*nCond+ic]);
				}
			}
		}
		if(((MaxAbs(nCond,satErr)<=etol)&&(MaxAbs(nCond,satCovErr)<=vtol))||(iter==maxIter))
		{
		  delete[] satCovErr;
		  if(nObv>0){
		    delete[] pointPtr;
		    delete bTrackPtr;
		    LaunchParticles(NULL,0,obvSet,nObv,iter);
		    CalcSaturationMoments(obvSet,0,nObv,oof,iter,0);
		  }
		  DeallocateIter();
		  break;
		}
		delete[] satCovErr;
		if (meas_err_opt == 1){
			for (int ierr1=0;ierr1<nCond;ierr1++){
				for (int ierr2=0;ierr2<nCond;ierr2++){
					satCovEP[ierr1*nCond+ierr2] = satCovEP[ierr1*nCond+ierr2]+swMeasErr[ierr1*nCond+ierr2];
				}
			}
		}
		//Calculate covariance of saturation and log-permeability
		CalcCovSYVector(iter,nCond,0);
		//Calculate conditioning weight (Lambda)
		CalcLambdaVector(nCond);
		//Update permeability moments
		UpdatePermeabilityMoments(0,nCond,iter);
		//Resolve flow problem
		flowPtr->deleteVeloPtr(); // To prevent memory leak by calling flowPtr->solve repeatedly
		flowPtr->solve(directory,is_flow);
		flowPtr->writeToFile(directory);
		domainPtr->deleteVeloCovariances(); // To prevent memory leak by calling flowPtr->writeToDomain repeatedly
		flowPtr->writeToDomain(domainPtr);
		//Dallocate HEAP
		DeallocateIter();
		iter++;
	}
	clockEnd = clock();
	delete krigPtr;
	cof<<endl<<"SATURATION CONDITIONING RUNTIME: "<<(clockEnd-clockStart)/CLOCKS_PER_SEC<<"  SECONDS";
	cof.close();
	oof.clear();
}

void SwCondiInverse::ResetCoordSetIndex(Coord* coordPtr)
{
	if(coordPtr != NULL)
	{
		coordPtr->index = NO_INDEX;
		ResetCoordSetIndex(coordPtr->coordPtr);
	}
}

int SwCondiInverse::CountMeasurements(int seqLoc)
{
	int nPoint = 1;
	for(int ii=seqLoc+1;ii<nCond;ii++)
	{
		if(condTime[ii]==condTime[seqLoc]){nPoint += 1;}
		else{break;}
	}
	return nPoint;
}

void SwCondiInverse::LaunchParticles(Coord** condPtr, int numCond, Coord** obvPtr, int numObv, int iter)
{
	ResetCoordSetIndex(coordSet);
	pointPtr = new Point[numCond+numObv];
	int nPoint = 0;
	for(int icond=0;icond<numCond;icond++)
	{
		if(condPtr[icond]->index==NO_INDEX)
		{
			condPtr[icond]->index = nPoint;
			Point tempPtr(condPtr[icond]->x_,
						  condPtr[icond]->y_,
						  0, //To be consistent
						  condPtr[icond]->i_,
						  condPtr[icond]->j_,
						  condPtr[icond]->k_);
			pointPtr[nPoint] = tempPtr;
			nPoint++;
		}
		condSwStrl[icond] = condPtr[icond]->index;
	}
	for(int iobv=0;iobv<numObv;iobv++)
	{
		if(obvPtr[iobv]->index==NO_INDEX)
		{
			obvPtr[iobv]->index = nPoint;
			Point tempPtr(obvPtr[iobv]->x_,
						  obvPtr[iobv]->y_,
						  0, //To be consistent
						  obvPtr[iobv]->i_,
						  obvPtr[iobv]->j_,
						  obvPtr[iobv]->k_);
			pointPtr[nPoint] = tempPtr;
			nPoint++;
		}
		obvSwStrl[iobv] = obvPtr[iobv]->index;
	}
	bTrackPtr = new ParticleTrack(domainPtr,nPoint,pointPtr);
	bTrackPtr->goTrackBackward();
	bTrackPtr->calcTrvlTimeMomentsNoWell();
	bTrackPtr->printTravelTime();
	/*// Display particle info
	char tauFileName[100];
	sprintf(tauFileName,"TravelTime_Iter=%d[satinv].out",iter);
	ofstream tof(tauFileName, ios::out);
	tof<<"xp yp distance tau_avg0tran tau_avg2tran tau_vartran tau_avg0carte tau_avg2carte tau_varcarte  tau_avg0 tau_avg2 tau_var"<<endl<<endl;
	for(int jj=0;jj<nPoint;jj++)
	{
		bTrackPtr->getParticle(jj)->PrintPathline(tof);
	}
	tof.close();*/
	//cout<<"Time "<<condTime[seqLoc]<<endl;
	//DisplayParticleInfo(numSw,bTrackPtr);
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

// condPointFlag indicates whether the data is conditioning or observing
void SwCondiInverse::CalcSaturationMoments(Coord** cPtrPtr, int pLoc, int numP, ofstream &pof, int iter, int condPointFlag)
{
	//Compute the moments
	Saturation* satPtr = new Saturation(bTrackPtr, viscosR, Swc, Sor, nS );
	if(condPointFlag==1){condSatPtr = satPtr;}
	else{obvSatPtr = satPtr;}
	//int obstInd;								//For Sequential Observation
	//(tInd==0) ? obstInd=1 : obstInd=0;			//For Sequential Observation
	//cout<<obstInd<<endl;						//For Sequential Observation
	//satPtr->solve(condiTimePtr[obstInd].ctime); //For Sequential Observation
	//satPtr->output((tInd+1)*10000+iter);		//For Sequential Observation
	//satPtr->solve(condTime[seqLoc]);
	int* swStrl;
	double* timeArray;
	if(condPointFlag==1)
	{
		swStrl = condSwStrl;
		timeArray = &condTime[pLoc];
	}
	else
	{
		swStrl = obvSwStrl;
		timeArray = &obvTime[pLoc];
	}
	int epf = 1;
	if(condPointFlag==0){epf = 0;}
	satPtr->solveMT(timeArray,numP,swStrl,epf);
	if(condPointFlag==1)
	{
		//satPtr->output((pLoc+1)*1000+iter); // Not necessary for conditioning
		//satCovEP = satPtr->calcSatCovEndPoints();

		//Extract saturation mean
		satVarEP = new double[numP];
		satAvgEP = new double[numP];
		satMeas = new double[numP];
		satErr = new double[numP];
	}
	if(condPointFlag==0)
	{satPtr->output(10000+iter);}
	for(int ind=0;ind<numP;ind++)
	{
		int sind = (ind*(((1-satPtr->endPointFlag)*satPtr->nNd)+satPtr->endPointFlag))+((1-satPtr->endPointFlag)*(bTrackPtr->getParticle(swStrl[ind])->getLength()-1));
		if(condPointFlag==1)
		{
			satVarEP[ind] = satPtr->getSatVar(sind);
			satAvgEP[ind] = satPtr->getSatAvg(sind);
			satMeas[ind] = swMeas[pLoc+ind];
			satErr[ind] = satMeas[ind]-satAvgEP[ind];
		}

		// Output to file
		pof<<"["<<timeArray[ind]<<","<<iter<<","<<ind<<"] ";
		if(condPointFlag==1){pof<<satMeas[ind]<<"  ";}
		pof<<satPtr->getSatAvg(sind)<<"  ";
		if(condPointFlag==1){pof<<satErr[ind]<<"  ";}
		pof<<satPtr->getSatVar(sind)<<endl;
	}
	pof<<endl;

	if(condPointFlag==1)
	{
		//Display mean
		cout<<endl<<"Sw Mean: ";
		for(int sind=0;sind<numP;sind++)
		{
			cout<<satAvgEP[sind]<<" ";
		}
		cout<<endl<<endl;
	}
}

void SwCondiInverse::DeallocateIter()
{
	delete[] condSwStrl;
	if(obvSwStrl!=NULL)
	  {
	    delete[] obvSwStrl;
	    obvSwStrl = NULL;
	  }
	delete[] pointPtr;
	delete bTrackPtr;
	delete condSatPtr;
	if(obvSatPtr!=NULL)
	  {
	    delete obvSatPtr;
	    obvSatPtr = NULL;
	  }
	if(satCovEP!=NULL)
	{
		delete[] satCovEP;
		satCovEP = NULL;
	}
	delete[] satVarEP;
	delete[] satAvgEP;
	delete[] satMeas;
	delete[] satErr;
	if(covSY!=NULL)
	{
		delete[] covSY;
		covSY = NULL;
		delete[] lambda;
		lambda = NULL;
	}
}

void SwCondiInverse::CalcCovSYVector(int iter, int numSw, int pLoc)
{
	covSY = new double[numSw*gridPtr->getNumPermNode()];

	char csyFileName[1000];
	sprintf(csyFileName,"C_SY_iter=%d.out",iter);
	//sprintf(csyFileName,"C_SY.out");


	ofstream csy(csyFileName, ios::out);
	csy<<"##### CSY #####"<<endl;

	for(int sInd=0;sInd<numSw;sInd++)
	{
		CollectDataS(condSwStrl[sInd]);
		for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
		{
			for(int pInd=0;pInd<bTrackPtr->getParticle(condSwStrl[sInd])->getLength();pInd++)
			{
				CalcCov_YVxi_YVeta_YEta(condSwStrl[sInd],yInd,pInd);
			}
			int syInd = sInd + (yInd*numSw);
			covSY[syInd] = CalcCovSYQuad(sInd,yInd,pLoc);
			csy<<sInd<<"   "<<yInd<<"   "<<covSY[syInd]<<"   "<<covSY[syInd]/sqrt(satVarEP[sInd]*permPtr->getYVar(yInd))<<endl;
		}
		DeallocateDataS();
	}
	csy.close();
}

void SwCondiInverse::CollectDataS(int lInd)
{
	int nPoint = bTrackPtr->getParticle(lInd)->getLength();
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
		ds[ii] = bTrackPtr->getParticle(lInd)->getDs(ii);
		dvdeta[ii] = bTrackPtr->getParticle(lInd)->getDvDeta(ii);
		costh[ii] = bTrackPtr->getParticle(lInd)->getCosth(ii);
		sinth[ii] = bTrackPtr->getParticle(lInd)->getSinth(ii);
		FindEdge(lInd,ii);
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
void SwCondiInverse::FindEdge(int lInd, int pInd)
{
	double xzero = 1e-6;

	x_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getX();
	y_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getY();
	z_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getZ();
	i_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getI();
	j_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getJ();
	k_[pInd] = bTrackPtr->getParticle(lInd)->getTrajectory(pInd).getK();

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

void SwCondiInverse::CalcCov_YVxi_YVeta_YEta(int lInd, int yInd, int pInd)
{
	double covYVxu,covYVxd,covYVyu,covYVyd; // u = upstream, d = downstream
	double covYVxlu,covYVxld,covYVxru,covYVxrd,covYVybu,covYVybd,covYVytu,covYVytd;
	// l = left, r = right, b = bottom, t = top
	RefreshPressureGradient();

	switch(edge[pInd])
	{
		case 0:
			covYVxlu = covYVxld = covYVxru = covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
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
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
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
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
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
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
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
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
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
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
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
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
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
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
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
				else if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]>=y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd]-1,j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>=x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]>=x_[pInd])&&(y_[pInd-1]>=y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<=x_[pInd+1])&&(y_[pInd]<=y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<=x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd]-1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]<=y_[pInd+1]))
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
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 5:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<x_[pInd])&&(y_[pInd-1]<=y_[pInd]))
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
				else if((x_[pInd-1]>=x_[pInd])&&(y_[pInd-1]<=y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]>=x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<=x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd]-1,j_[pInd]+1,k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]<=x_[pInd+1])&&(y_[pInd]>=y_[pInd+1]))
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
				else if((x_[pInd]>x_[pInd+1])&&(y_[pInd]>=y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]-1,j_[pInd],k_[pInd],i_[pInd]-1,j_[pInd]+1,k_[pInd],yInd);
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 6:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<=x_[pInd])&&(y_[pInd-1]<y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else if((x_[pInd-1]<=x_[pInd])&&(y_[pInd-1]>=y_[pInd]))
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
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]>=y_[pInd]))
				{
					covYVxlu = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxru = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybu = CalcCovYVy(i_[pInd]+1,j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
			{covYVxld = covYVxrd = covYVybd = covYVytd = 0;}
			else
			{
				if((x_[pInd]<x_[pInd+1])&&(y_[pInd]<=y_[pInd+1]))
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
				else if((x_[pInd]>=x_[pInd+1])&&(y_[pInd]<=y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>=x_[pInd+1])&&(y_[pInd]>y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd]+1,j_[pInd]-1,k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd]-1,k_[pInd],i_[pInd],j_[pInd],k_[pInd],yInd);
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			break;
		case 7:
			if(pInd == 0){covYVxlu = covYVxru = covYVybu = covYVytu = 0;}
			else
			{
				if((x_[pInd-1]<=x_[pInd])&&(y_[pInd-1]<=y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybu = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd-1]<=x_[pInd])&&(y_[pInd-1]>y_[pInd]))
				{
					covYVxlu = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxru = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVybu = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytu = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd-1]>x_[pInd])&&(y_[pInd-1]<=y_[pInd]))
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
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
			}
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
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
				else if((x_[pInd]<x_[pInd+1])&&(y_[pInd]>=y_[pInd+1]))
				{
					covYVxld = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVxrd = 0; // Won't be use in interpulation (the point is at a corner on left edge)
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd]+1,j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
				}
				else if((x_[pInd]>=x_[pInd+1])&&(y_[pInd]<y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd]+1,k_[pInd],i_[pInd]+1,j_[pInd]+1,k_[pInd],yInd);
					covYVybd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
					covYVytd = 0; // Won't be use in interpulation (the point is at a corner on bottom edge)
				}
				else if((x_[pInd]>=x_[pInd+1])&&(y_[pInd]>=y_[pInd+1]))
				{
					covYVxld = 0; // Won't be use in interpulation (the point is at a corner on right edge)
					covYVxrd = CalcCovYVx(i_[pInd],j_[pInd],k_[pInd],i_[pInd]+1,j_[pInd],k_[pInd],yInd);
					covYVybd = 0; // Won't be use in interpulation (the point is at a corner on top edge)
					covYVytd = CalcCovYVy(i_[pInd],j_[pInd],k_[pInd],i_[pInd],j_[pInd]+1,k_[pInd],yInd);
				}
				else{cerr<<"Error in directions, check SwCondiInverse.cc"<<endl;}
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
			if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
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

	//if((covYVxlu!=0)&&(covYVxru!=0))
	//{
		covYVxu = ((covYVxlu*(xrt_[pInd]-x_[pInd]))+(covYVxru*(x_[pInd]-xlb_[pInd])))/(xrt_[pInd]-xlb_[pInd]);
	//}else{covYVxu = 0;}
	//if((covYVxld!=0)&&(covYVxrd!=0))
	//{
		covYVxd = ((covYVxld*(xrt_[pInd]-x_[pInd]))+(covYVxrd*(x_[pInd]-xlb_[pInd])))/(xrt_[pInd]-xlb_[pInd]);
	//}else{covYVxd = 0;}
	//if((covYVybu!=0)&&(covYVytu!=0))
	//{
		covYVyu = ((covYVybu*(yrt_[pInd]-y_[pInd]))+(covYVytu*(y_[pInd]-ylb_[pInd])))/(yrt_[pInd]-ylb_[pInd]);
	//}else{covYVyu = 0;}
	//if((covYVybd!=0)&&(covYVytd!=0))
	//{
		covYVyd = ((covYVybd*(yrt_[pInd]-y_[pInd]))+(covYVytd*(y_[pInd]-ylb_[pInd])))/(yrt_[pInd]-ylb_[pInd]);
	//}else{covYVyd = 0;}

	if(pInd == 0){covYVxiU[pInd] = covYVetaU[pInd] = 0;}
	else
	{
		covYVxiU[pInd] = (covYVxu*costh[pInd])+(covYVyu*sinth[pInd]);
		covYVetaU[pInd] = -(covYVxu*sinth[pInd])+(covYVyu*costh[pInd]);
	}
	if(pInd == bTrackPtr->getParticle(lInd)->getLength()-1)
	{covYVxiD[pInd] = covYVetaD[pInd] = 0;}
	else
	{
		covYVxiD[pInd] = (covYVxd*costh[pInd])+(covYVyd*sinth[pInd]);
		covYVetaD[pInd] = -(covYVxd*sinth[pInd])+(covYVyd*costh[pInd]);
	}
	if(pInd == 0){covYEta[pInd] = 0;}
	else
	{
		double vxiavg_m1 = bTrackPtr->getParticle(lInd)->getVAvg(pInd-1);
		double vxiavg = bTrackPtr->getParticle(lInd)->getVAvg(pInd);
		covYEta[pInd] = covYEta[pInd-1] + (covYVetaD[pInd-1]/vxiavg_m1*ds[pInd]/2) + (covYVetaU[pInd]/vxiavg*ds[pInd]/2);
	}
}

void SwCondiInverse::RefreshPressureGradient()
{
	dpdx = P0e->getDP0DXi(0,0);
	dpdy = P0e->getDP0DXi(1,0);
}

double SwCondiInverse::CalcCovYVx(int i1, int j1, int k1, int i2, int j2, int k2, int yInd)
{//To have proper input, i1 = i2-1, j1 = j2, k1 = k2,
	double perm_interface[3];
	int permInd, ijkpermtemp_;
	double dpdxtemp, cyytemp, dcpydxtemp;
	int ijk1, ijk2;
	double poro;

	dpdxtemp = dcpydxtemp = 0; // Will not be changed if the interface is at a border (no flow)
	ijk1 = gridPtr->getIndex(i1,j1,k1);
	ijk2 = gridPtr->getIndex(i2,j2,k2);

	if(i1>=0)
	{
		poro = domainPtr->getPoro(i1,j1,k1);
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
		poro = domainPtr->getPoro(i2,j2,k2);
		permPtr->Perm_x(i2,j2,k2,perm_interface);
		permInd = 0;
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i2)-1,gridPtr->toJPerm(j2),gridPtr->toKPerm(k2));
	}

	cyytemp = permPtr->getCYY(ijkpermtemp_,yInd);
	return(-perm_interface[permInd]/poro*((dpdxtemp*cyytemp)+dcpydxtemp));
}

double SwCondiInverse::CalcCovYVy(int i1, int j1, int k1, int i2, int j2, int k2, int yInd)
{//To have proper input, i1 = i2, j1 = j2-1, k1 = k2
	double perm_interface[3];
	int permInd, ijkpermtemp_;
	double dpdytemp, cyytemp, dcpydytemp;
	int ijk1, ijk2;
	double poro = 1;

	dpdytemp = dcpydytemp = 0; // Will not be changed if the interface is at a border (no flow)
	ijk1 = gridPtr->getIndex(i1,j1,k1);
	ijk2 = gridPtr->getIndex(i2,j2,k2);

	if(j1>=0)
	{
		poro = domainPtr->getPoro(i1,j1,k1);
		permPtr->Perm_y(i1,j1,k1,perm_interface);
		permInd = 2;
		if(j1<(gridPtr->getNy()-1))
		{
			dpdytemp = dpdy[gridPtr->getIndex(i1,j1,k1)];
			dcpydytemp = (CYPe->getCYP(ijk2,yInd)-CYPe->getCYP(ijk1,yInd))/gridPtr->getDy(j1);
		}
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i1),gridPtr->toJPerm(j1)+1,gridPtr->toKPerm(k1));
	}
	else
	{
		poro = domainPtr->getPoro(i2,j2,k2);
		permPtr->Perm_y(i2,j2,k2,perm_interface);
		permInd = 0;
		ijkpermtemp_ = gridPtr->getPermIndex(gridPtr->toIPerm(i2),gridPtr->toJPerm(j2)-1,gridPtr->toKPerm(k2));
	}

	cyytemp = permPtr->getCYY(ijkpermtemp_,yInd);
	return(-perm_interface[permInd]/poro*((dpdytemp*cyytemp)+dcpydytemp));
}

double SwCondiInverse::CalcCovSYQuad(int sInd, int yInd, int pLoc)
{
	double savg = satAvgEP[sInd];
	double yavg = permPtr->getYAvg(yInd);
	double yvar = permPtr->getYVar(yInd);
	double tauavg = bTrackPtr->getParticle(condSwStrl[sInd])->getLastTravelTime();
	double tauvar = bTrackPtr->getParticle(condSwStrl[sInd])->getLastTravelTVar();
	double lntauavg, lntauvar;
	int lenn = 1;
	transform1(lenn,&tauavg,&tauvar,&lntauavg,&lntauvar);
	double covYTau = CalcCovYTau(condSwStrl[sInd]);
	double covYLnTau = CalcCovYLnTau(tauavg,covYTau);
	double rhoYLnTau = covYLnTau/sqrt(yvar*lntauvar);

	//2D quadrature integration
	double sDRange = 3.5, intRange = 1;
	int qPoint = 3;
	double* rloc, * weight;
	GaussianQuadrature(qPoint,rloc,weight);

	double ymin = yavg - (sDRange*sqrt(yvar));
	double ymax = yavg + (sDRange*sqrt(yvar));
	double lntaumin = lntauavg - (sDRange*sqrt(lntauvar));
	double lntaumax = lntauavg + (sDRange*sqrt(lntauvar));

	double covSY = 0;
	double lPointY,rPointY,lPointLntau,rPointLntau;
	double lntauFront = log(condSatPtr->fpwStar*condTime[pLoc+sInd]);

	lPointY = ymin;
	while(true)
	{
		rPointY = lPointY+(intRange*sqrt(yvar));
		if(rPointY>ymax){rPointY = ymax;}
		for(int iy=0;iy<qPoint;iy++)
		{
			double yVal = ((rPointY-lPointY)/2*rloc[iy])+((rPointY+lPointY)/2);
			lPointLntau = lntaumin;
			while(true)
			{
				rPointLntau = lPointLntau+(intRange*sqrt(lntauvar));
				if(rPointLntau>lntaumax){rPointLntau = lntaumax;}
				if((lntauFront<rPointLntau)&&(lntauFront>lPointLntau)){rPointLntau = lntauFront;}
				for(int ilt=0;ilt<qPoint;ilt++)
				{
					double lntauVal = ((rPointLntau-lPointLntau)/2*rloc[ilt])+((rPointLntau+lPointLntau)/2);
					double tauVal = exp(lntauVal);
					double sVal;
					condSatPtr->t = condTime[pLoc+sInd];
					condSatPtr->calcSatDet(1,&tauVal,&sVal);
					double lntsd = (lntauVal-lntauavg)/sqrt(lntauvar);
					double ysd = (yVal-yavg)/sqrt(yvar);
					double param = (sVal-savg)*(yVal-yavg);
					covSY += ((rPointY-lPointY)/2)*weight[iy]*((rPointLntau-lPointLntau)/2)*weight[ilt]*
						     param*(1/(2*3.1415926535897932*sqrt(lntauvar*yvar*(1-(rhoYLnTau*rhoYLnTau))))*
							 exp(-0.5/(1-(rhoYLnTau*rhoYLnTau))*((lntsd*lntsd)-(2*rhoYLnTau*lntsd*ysd)+(ysd*ysd))));
				}
				if(rPointLntau==lntaumax){break;}
				lPointLntau = rPointLntau;
			}
		}
		if(rPointY==ymax){break;}
		lPointY = rPointY;
	}
	delete[] rloc;
	delete[] weight;
	return covSY;
}

double SwCondiInverse::CalcCovYTau(int lInd)
{
	double covYTau = 0;
	for(int pInd=1;pInd<bTrackPtr->getParticle(lInd)->getLength();pInd++)
	{
		double vxiavg = bTrackPtr->getParticle(lInd)->getVAvg(pInd-1);
		//Use dvdeta[pInd] instead of dvdeta[pInd-1] because dvdeta is cell-center variable
		covYTau += 1/(vxiavg*vxiavg)*(covYVxiD[pInd-1]+(covYEta[pInd-1]*dvdeta[pInd]))*ds[pInd]/2;
		vxiavg = bTrackPtr->getParticle(lInd)->getVAvg(pInd);
		covYTau += 1/(vxiavg*vxiavg)*(covYVxiU[pInd]+(covYEta[pInd]*dvdeta[pInd]))*ds[pInd]/2;
	}
	return -covYTau;
}

double SwCondiInverse::CalcCovYLnTau(double tauMean, double covYTau)
{
	//This implementation is to test analytical convertion, replacing the above computation
	//double lntauMean = (2*log(tauMean))-(log(pow(tauMean,2)+tauVar)/2);
	//double lntauVar = log(1+(tauVar/pow(tauMean,2)));
	//return((1/exp((lntauVar/2)+lntauMean)*(covYTau+(yMean*tauMean)))-yMean);
	return(covYTau/tauMean);
}

void SwCondiInverse::CalcLambdaVector(int numSw)
{
	lambda = new double[numSw*gridPtr->getNumPermNode()];
	double* satCovEPInv = new double[numSw*numSw];
	krigPtr->calCovInverse(numSw,satCovEP,satCovEPInv);

	for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
	{
		for(int sInd1=0;sInd1<numSw;sInd1++)
		{
			lambda[sInd1+(yInd*numSw)] = 0;
			for(int sInd2=0;sInd2<numSw;sInd2++)
			{
				lambda[sInd1+(yInd*numSw)] += satCovEPInv[sInd2+(sInd1*numSw)]*covSY[sInd2+(yInd*numSw)];
			}
			lambda[sInd1+(yInd*numSw)] = dFactor*lambda[sInd1+(yInd*numSw)];
		}
	}
	delete[] satCovEPInv; // Used to leak
}

void SwCondiInverse::UpdatePermeabilityMoments(int seqLoc, int numSw, int iter)
{
	int numY = gridPtr->getNumPermNode();
	char ymFileName[1000], yvFileName[1000];
	sprintf(ymFileName,"YMean_Time=%.3f_Iter=%d.out",condTime[seqLoc],iter+1);
	ofstream ym(ymFileName, ios::out);
	ym<<"##### Y Mean #####"<<endl;
	sprintf(yvFileName,"YVar_Time=%.3f_Iter=%d.out",condTime[seqLoc],iter+1);
	ofstream yv(yvFileName, ios::out);
	yv<<"##### Y Variance #####"<<endl;

	for(int yInd1=0;yInd1<numY;yInd1++)
	{
		//Mean update
		double meann = permPtr->getYAvg(yInd1);
		for(int sInd=0;sInd<numSw;sInd++)
		{
			meann += lambda[sInd+(yInd1*numSw)]*satErr[sInd];
		}
		permPtr->setYAvg(yInd1,meann);
		ym<<yInd1<<"   "<<meann<<endl;


		//CYY update
		for(int yInd2=yInd1;yInd2<numY;yInd2++)
		{
			double lambdaCSY = 0;
			for(int sInd=0;sInd<numSw;sInd++)
			{
				lambdaCSY += lambda[sInd+(yInd1*numSw)]*covSY[sInd+(yInd2*numSw)];
				if (meas_err_opt == 2){
					for(int sInd2=0;sInd2<numSw;sInd2++){
						lambdaCSY -= lambda[sInd+(yInd1*numSw)]*swMeasErr[(sInd+seqLoc)*nCond+(sInd2+seqLoc)]*lambda[sInd2+(yInd2*numSw)];
					}
				}
			}
			permPtr->setCYY_Cond(yInd1,yInd2,lambdaCSY);
			if(yInd1!=yInd2){permPtr->setCYY_Cond(yInd2,yInd1,lambdaCSY);}
			//else{yv<<yInd1<<"   "<<permPtr->getCYY(yInd1,yInd2)-lambdaCSY<<endl;} //WHY -lambdaCSY????
			else{yv<<yInd1<<"   "<<permPtr->getCYY(yInd1,yInd2)<<endl;}
		}
	}
	//ofstream covyall("CYY_ALL.out",ios::out);
	//ofstream coryall("CorrYY_ALL.out",ios::out);
	/*for(int yInd1=0;yInd1<numY;yInd1++)
	{
		for(int yInd2=yInd1;yInd2<numY;yInd2++)
		{
			double corryy = permPtr->getCYY(yInd1,yInd2)/sqrt(permPtr->getCYY(yInd1,yInd1)*permPtr->getCYY(yInd2,yInd2));
			//covyall<<permPtr->getCYY(yInd1,yInd2)<<" ";
			//coryall<<corryy<<" ";
			if(sqrt(corryy*corryy)>1)
			{
				cout<<"This is WRONG corryy cannot greater than 1"<<endl;
			}
		}
		//covyall<<endl;
		//coryall<<endl;
	}
	//covyall.close();
	//coryall.close();*/
	ym.close();
	yv.close();
	permPtr->wrtPermMomnt();
}

/*double SwCondiInverse::CalcCovSY(int sInd, int yInd)
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
	double covYLnTau = CalcCovYLnTau(yavg,tauavg,tauvar,covYTau);
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
}*/
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

