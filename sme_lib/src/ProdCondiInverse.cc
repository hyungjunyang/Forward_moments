#include "ProdCondiInverse.h"

ProdCondiInverse::ProdCondiInverse(Domain* domain, Flow* flow)
: domainPtr(domain), flowPtr(flow)
{
	wellsPtr = domainPtr->getWellArr();
	gridPtr = domainPtr->getGrid();
	permPtr = flowPtr->getPermPtr();
	Swc     = domainPtr->getSwc();
	Sor     = domainPtr->getSor();
	viscosR = domainPtr->getViscosR();
	blSolnPtr = new BLSolution(1000, Swc,  Sor, viscosR);
	prodPtr = new Production(domainPtr->getUnit(), domainPtr, flowPtr, blSolnPtr, 1, 1, false);
	swCondPtr = new SwCondiInverse(domainPtr,flowPtr);

	//coordSet = NULL;
	condSet = NULL;
	obvSet = NULL;
	prodCov = NULL;
	prodVar = NULL;
	prodAvg = NULL;
	prodErr = NULL;
	covQY = NULL;
	lambda = NULL;

	nMeas = nObv = 0;
	readDataFlag = 0;

	cout<<"ProdCondiInverse object created."<<endl;
}

ProdCondiInverse::~ProdCondiInverse()
{
	delete blSolnPtr;
	delete prodPtr;
	delete swCondPtr;
	if(readDataFlag==1)
	{
		if(nWCond>0)
		{
			delete[] condSet;
			delete[] condQtMeas;
			for(int ii=0;ii<nWCond;ii++)
			{
				if(condQwTime[ii]!=NULL)
				{
					delete[] condQwTime[ii];
					delete[] condQwMeas[ii];
				}
			}
			delete[] condQwTime;
			delete[] condQwMeas;
		}
		if(nWObv>0)
		{
			delete[] obvSet;
			for(int ii=0;ii<nWObv;ii++)
			{
				if(obvQwTime[ii]!=NULL)
				{
					delete[] obvQwTime[ii];
				}
			}
			delete[] obvQwTime;
		}
	}
}

void ProdCondiInverse::ProductionConditioning(char *directory, ifstream &is_flow)
{
	clockStart = clock();
	ReadData();
	prodPtr->SetParticlesAtProdWells();
	cout<<"Start Production Conditioning"<<endl;

	krigPtr = new KrigWeigh(gridPtr,permPtr,nMeas);
	// Output files
	char prodFileName[100];
	sprintf(prodFileName,"CondProdMoments[prodinv].out");
	ofstream cof(prodFileName, ios::out);
	cof<<"[Well#,Iteration#,Type] Time Measurement Prediction Diff Variance"<<endl<<endl;
	sprintf(prodFileName,"ObvProdMoments[prodinv].out");
	ofstream oof(prodFileName, ios::out);
	oof<<"[Well#,Iteration#,Type] Time Prediction Variance"<<endl<<endl;

	//Iteration
	int iter = 0;
	while(true)
	{
		cout<<"Iteration "<<iter<<endl;

		//Calculate Production Moments
		prodPtr->LaunchParticles();
		cout<<"Done - LaunchParticles"<<endl;
		prodPtr->ClearCovYLnTau();
		cout<<"Done - ClearCovYLnTau()"<<endl;
		if(nWCond>0){CalcProductionMoments(cof,iter,1);}
		cout<<"Done - CalcProductionMoments for Cond"<<endl;
		//exit(0);

		//Calculate Observation at the end only
		if(((MaxAbs(nMeas,prodErr)<=etol)&&(MaxAbs(nMeas,prodVar)<=vtol))||(iter==maxIter)||(iter==0))
		{
		  if(nWObv>0){CalcProductionMoments(oof,iter,0);}
		  cout<<"Done - CalcProductionMoments for Obv"<<endl;
		}

		//Check Convergence Criteria
		if(((MaxAbs(nMeas,prodErr)<=etol)&&(MaxAbs(nMeas,prodVar)<=vtol))||(iter==maxIter))
		{
			DeallocateIter();
			break;
		}

		//Calculate LHS for Kriging
		CalcCovQQMatrix();
		cout<<"Done - CalcCovQQMatrix()"<<endl;

		//Calculate RHS for Kriging
		CalcCovQYVector();
		//break;
		cout<<"Done - CalcCovQYVector()"<<endl;

		//Calculate Weight Vector (Lambda)
		CalcLambdaVector();
		cout<<"Done - CalcLmabdaVector"<<endl;

		//Update Permeability Field
		UpdatePermeabilityMoments(iter);
		cout<<"Done - UdatePermeabilityMoments()"<<endl;

		//Resolve for Flow Problem
		flowPtr->deleteVeloPtr(); // To prevent memory leak by calling flowPtr->solve repeatedly
		cout<<"Start - solving flow problem"<<endl;
		flowPtr->solve(directory,is_flow);
		cout<<"Done -  solving flow problem"<<endl;
		flowPtr->writeToFile(directory);
		domainPtr->deleteVeloCovariances(); // To prevent memory leak by calling flowPtr->writeToDomain repeatedly
		flowPtr->writeToDomain(domainPtr);

		//Deallocate HEAP
		DeallocateIter();
		cout<<"Done - deallocating"<<endl;

		iter++;
	}
	delete krigPtr;
	clockEnd = clock();
	cof<<endl<<"PRODUCTION CONDITIONING RUNTIME: "<<(clockEnd-clockStart)/CLOCKS_PER_SEC<<"  SECONDS";
	cof.close();
	oof.close();
}

void ProdCondiInverse::ReadData()
{
	cout<<"Start reading data for ProdCondiInverse"<<endl;
	ifstream in_ProdCondi;
	in_ProdCondi.open("prod_condi.in",ios::in);
	if(in_ProdCondi.bad())
	{
		cerr<<"Can not open file, prod_condi.in"<<endl;
		exit(8);
	}
	else
	{
		cerr<<"File prod_condi.in was opened successfully!"<<endl;
	}

	Junk(in_ProdCondi);
	in_ProdCondi >> dFactor;
	Junk(in_ProdCondi);
	in_ProdCondi >> nWCond >> nWObv;
	if(nWCond>0)
	{
		condSet = new WCondData[nWCond];
		condQtMeas = new double[nWCond];
		condQwTime = new double*[nWCond];
		condQwMeas = new double*[nWCond];
		for(int ii=0;ii<nWCond;ii++)
		{
			Junk(in_ProdCondi);
			in_ProdCondi >> condSet[ii].wIndex >> condSet[ii].calcQtFlag >> condSet[ii].numQw;
			condSet[ii].wIndex--; //Back to "0" based index
			if(condSet[ii].wIndex>=domainPtr->getNumWells())
			{
				cerr<<"Well "<<condSet[ii].wIndex+1<<" does not exist! There are total of "
					<<domainPtr->getNumWells()<<" wells";
				exit(8);
			}
			if(condSet[ii].calcQtFlag==1)
			{
				nMeas++;
				Junk(in_ProdCondi);
				in_ProdCondi >> condQtMeas[ii];
			}
			if(condSet[ii].numQw>0)
			{
				nMeas += condSet[ii].numQw;
				prodPtr->ProdWellIs(condSet[ii].wIndex);
				condQwTime[ii] = new double[condSet[ii].numQw];
				condQwMeas[ii] = new double[condSet[ii].numQw];
				for(int jj=0;jj<condSet[ii].numQw;jj++)
				{
					Junk(in_ProdCondi);
					in_ProdCondi >> condQwTime[ii][jj] >> condQwMeas[ii][jj];
				}
			}
			else
			{
				condQwTime[ii] = NULL;
				condQwMeas[ii] = NULL;
			}
		}
	}
	if(nWObv>0)
	{
		obvSet = new WCondData[nWObv];
		obvQwTime = new double*[nWObv];
		for(int ii=0;ii<nWObv;ii++)
		{
			Junk(in_ProdCondi);
			in_ProdCondi >> obvSet[ii].wIndex >> obvSet[ii].calcQtFlag >> obvSet[ii].numQw;
			obvSet[ii].wIndex--; // Back to "0" based index
			nObv += obvSet[ii].calcQtFlag;
			if(obvSet[ii].numQw>0)
			{
				nObv += obvSet[ii].numQw;
				prodPtr->ProdWellIs(obvSet[ii].wIndex);
				obvQwTime[ii] = new double[obvSet[ii].numQw];
				for(int jj=0;jj<obvSet[ii].numQw;jj++)
				{
					Junk(in_ProdCondi);
					in_ProdCondi >> obvQwTime[ii][jj];
				}
			}
			else
			{
				obvQwTime[ii] = NULL;
			}
		}
	}
	Junk(in_ProdCondi);
	in_ProdCondi >> maxIter >> etol >> vtol;
	in_ProdCondi.close();
	readDataFlag = 1;
}

void ProdCondiInverse::CalcProductionMoments(ofstream &pof, int iter, int condWellFlag)
{
	int iMeas = 0;
	int nWell;
	WCondData* calcSet;
	double** calcQwTime;
	double* calcAvg;
	double* calcVar;
	if(condWellFlag==1)
	{
		nWell = nWCond;
		calcSet = condSet;
		calcQwTime = condQwTime;
		prodAvg = new double[nMeas];
		prodErr = new double[nMeas];
		prodVar = new double[nMeas];
		calcAvg = prodAvg;
		calcVar = prodVar;
	}
	else
	{
		nWell = nWObv;
		calcSet = obvSet;
		calcQwTime = obvQwTime;
		calcAvg = new double[nObv];
		calcVar = new double[nObv];
	}
	for(int iw=0;iw<nWell;iw++)
	{
		if(calcSet[iw].calcQtFlag==1)
		{
			pof<<"["<<calcSet[iw].wIndex+1<<","<<iter<<",Qt] N/A  ";
			calcAvg[iMeas] = prodPtr->CalcTotalProdRateMean(calcSet[iw].wIndex);
			if(condWellFlag==1)
			{
				prodErr[iMeas] = condQtMeas[iw]-calcAvg[iMeas];
				pof<<condQtMeas[iw]<<"  ";
			}
			pof<<calcAvg[iMeas]<<"  ";
			if(condWellFlag==1){pof<<prodErr[iMeas]<<"  ";}
			calcVar[iMeas] = prodPtr->CalcTotalProdRateCov(calcSet[iw].wIndex,calcSet[iw].wIndex);
			pof<<calcVar[iMeas]<<endl;
			iMeas++;
		}
		if(calcSet[iw].numQw>0)
		{
			for(int ii=0;ii<calcSet[iw].numQw;ii++)
			{
				pof<<"["<<calcSet[iw].wIndex+1<<","<<iter<<",Qw] "<<calcQwTime[iw][ii]<<"  ";
				calcAvg[iMeas] = prodPtr->CalcWaterProdRateMean(calcSet[iw].wIndex,calcQwTime[iw][ii]);
				if(condWellFlag)
				{
					prodErr[iMeas] = condQwMeas[iw][ii]-calcAvg[iMeas];
					pof<<condQwMeas[iw][ii]<<"  ";
				}
				pof<<calcAvg[iMeas]<<"  ";
				if(condWellFlag==1){pof<<prodErr[iMeas]<<"  ";}
				calcVar[iMeas] = prodPtr->CalcWaterProdRateCov(calcSet[iw].wIndex,calcQwTime[iw][ii],calcAvg[iMeas],calcSet[iw].wIndex,calcQwTime[iw][ii],calcAvg[iMeas]);
				//double qtqw = prodPtr->CalcTotal_WaterProdRateCov(calcSet[iw].wIndex,calcSet[iw].wIndex,calcQwTime[iw][ii]);
				//pof<<calcVar[iMeas]<<" "<<qtqw<<endl;
				pof<<calcVar[iMeas]<<endl;
				iMeas++;
			}
		}
	}
	pof<<endl;
	if(condWellFlag==0)
	{
		delete[] calcAvg;
		delete[] calcVar;
	}
}

void ProdCondiInverse::DeallocateIter()
{
	delete[] prodAvg;
	delete[] prodErr;
	delete[] prodVar;
	if(prodCov!=NULL)
	{
		delete[] prodCov;
		prodCov = NULL;
	}
	if(covQY!=NULL)
	{
		delete[] covQY;
		covQY = NULL;
	}
	if(lambda!=NULL)
	{
		delete[] lambda;
		lambda = NULL;
	}
}

void ProdCondiInverse::CalcCovQQMatrix()
{
	prodCov = new double[nMeas*nMeas];
	int iMeas1 = 0;
	for(int wCInd1=0;wCInd1<nWCond;wCInd1++)
	{
		int wInd1 = condSet[wCInd1].wIndex;
		if(condSet[wCInd1].calcQtFlag==1)
		{
			int iMeas2 = iMeas1;
			for(int wCInd2=wCInd1;wCInd2<nWCond;wCInd2++)
			{
				int wInd2 = condSet[wCInd2].wIndex;
				if(condSet[wCInd2].calcQtFlag==1)
				{
					if(iMeas2==iMeas1){prodCov[iMeas2+(iMeas1*nMeas)] = prodVar[iMeas1];}
					else
					{
						prodCov[iMeas2+(iMeas1*nMeas)] = prodPtr->CalcTotalProdRateCov(wInd1,wInd2);
						prodCov[iMeas1+(iMeas2*nMeas)] = prodCov[iMeas2+(iMeas1*nMeas)];
					}
					iMeas2++;
				}
				if(condSet[wCInd2].numQw>0)
				{
					for(int iqw2=0;iqw2<condSet[wCInd2].numQw;iqw2++)
					{
						prodCov[iMeas2+(iMeas1*nMeas)] = 
						  prodPtr->CalcTotal_WaterProdRateCov(wInd1,wInd2,condQwTime[wCInd2][iqw2]);
						prodCov[iMeas1+(iMeas2*nMeas)] = prodCov[iMeas2+(iMeas1*nMeas)];
						iMeas2++;
					}
				}
			}
			iMeas1++;
		}
		if(condSet[wCInd1].numQw>0)
		{
			for(int iqw1=0;iqw1<condSet[wCInd1].numQw;iqw1++)
			{
				int iMeas2 = iMeas1;
				for(int wCInd2=wCInd1;wCInd2<nWCond;wCInd2++)
				{
					int wInd2 = condSet[wCInd2].wIndex;
					if((condSet[wCInd2].calcQtFlag==1)&&(wCInd2>wCInd1))
					{
						prodCov[iMeas2+(iMeas1*nMeas)] = 
						  prodPtr->CalcTotal_WaterProdRateCov(wInd2,wInd1,condQwTime[wCInd1][iqw1]);
						prodCov[iMeas1+(iMeas2*nMeas)] = prodCov[iMeas2+(iMeas1*nMeas)];
						iMeas2++;
					}
					if(condSet[wCInd2].numQw>0)
					{
						int iqw2Start;
						if(wCInd2==wCInd1){iqw2Start = iqw1;}
						else{iqw2Start = 0;}
						for(int iqw2=iqw2Start;iqw2<condSet[wCInd2].numQw;iqw2++)
						{
							if(iMeas2==iMeas1){prodCov[iMeas2+(iMeas1*nMeas)] = prodVar[iMeas1];}
							else
							{
							  prodCov[iMeas2+(iMeas1*nMeas)] = 
							    prodPtr->CalcWaterProdRateCov(wInd1,condQwTime[wCInd1][iqw1],prodAvg[iMeas1],
											  wInd2,condQwTime[wCInd2][iqw2],prodAvg[iMeas2]);
							  prodCov[iMeas1+(iMeas2*nMeas)] = prodCov[iMeas2+(iMeas1*nMeas)];
							}
							iMeas2++;
						}
					}
				}
				iMeas1++;
			}
		}
	}
	cout<<endl<<"C_{QQ}"<<endl;
	for(int ii=0;ii<nMeas;ii++)
	{
		for(int jj=0;jj<nMeas;jj++)
		{
			cout<<prodCov[ii+(jj*nMeas)]<<"  ";
		}
		cout<<endl;
	}
}

void ProdCondiInverse::CalcCovQYVector()
{
	covQY = new double[nMeas*gridPtr->getNumPermNode()];
	ofstream cqy("C_QY.out", ios::out);
	cqy<<"##### CQY #####"<<endl;

	swCondPtr->RefreshPressureGradient();

	int iMeas = 0;
	for(int wCInd=0;wCInd<nWCond;wCInd++)
	{
		int wInd = condSet[wCInd].wIndex;
		int wLen = wellsPtr[wInd].getWellLen();
		if(condSet[wCInd].calcQtFlag==1)
		{
			for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
			{
				int qyInd = iMeas + (yInd*nMeas);
				covQY[qyInd] = 0;
				for(int ib=0;ib<wLen;ib++)
				{
					int ii = wellsPtr[wInd].getWellIJK(0,ib);
					int jj = wellsPtr[wInd].getWellIJK(1,ib);
					int kk = wellsPtr[wInd].getWellIJK(2,ib);
					double dx = domainPtr->getBlock(ii,jj,kk)->getDx();
					double dy = domainPtr->getBlock(ii,jj,kk)->getDy();
					double dz = domainPtr->getBlock(ii,jj,kk)->getDz();
					double poro = domainPtr->getPoro(ii,jj,kk);
					covQY[qyInd] += swCondPtr->CalcCovYVx(ii-1,jj,kk,ii,jj,kk,yInd)*poro*dy*dz;
					covQY[qyInd] -= swCondPtr->CalcCovYVx(ii,jj,kk,ii+1,jj,kk,yInd)*poro*dy*dz;
					covQY[qyInd] += swCondPtr->CalcCovYVy(ii,jj-1,kk,ii,jj,kk,yInd)*poro*dx*dz;
					covQY[qyInd] -= swCondPtr->CalcCovYVy(ii,jj,kk,ii,jj+1,kk,yInd)*poro*dx*dz;
				}
				cqy<<iMeas<<"   "<<yInd<<"   "<<covQY[qyInd]<<endl;
			}
			iMeas++;
		}
		if(condSet[wCInd].numQw>0)
		{
			for(int iqw=0;iqw<condSet[wCInd].numQw;iqw++)
			{
				for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
				{
				        int qyInd = iMeas + (yInd*nMeas);
				        covQY[qyInd] = prodPtr->CalcCrossCovQwY(condSet[wCInd].wIndex,condQwTime[wCInd][iqw],yInd);
					cqy<<iMeas<<"   "<<yInd<<"   "<<covQY[qyInd]<<endl;
				        /*
					for(int ib=0;ib<wLen;ib++)
					{
					}
					*/
				}
				iMeas++;
			}
		}
	}
	cqy.close();
}

void ProdCondiInverse::CalcLambdaVector()
{
	lambda = new double[nMeas*gridPtr->getNumPermNode()];
	double* prodCovInv = new double[nMeas*nMeas];
	krigPtr->calCovInverse(nMeas,prodCov,prodCovInv);

	for(int yInd=0;yInd<gridPtr->getNumPermNode();yInd++)
	{
		for(int qInd1=0;qInd1<nMeas;qInd1++)
		{
			lambda[qInd1+(yInd*nMeas)] = 0;
			for(int qInd2=0;qInd2<nMeas;qInd2++)
			{
				lambda[qInd1+(yInd*nMeas)] += prodCovInv[qInd2+(qInd1*nMeas)]*covQY[qInd2+(yInd*nMeas)];
			}
			lambda[qInd1+(yInd*nMeas)] = dFactor*lambda[qInd1+(yInd*nMeas)];
		}
	}
	delete[] prodCovInv;
}

void ProdCondiInverse::UpdatePermeabilityMoments(int iter)
{
	int numY = gridPtr->getNumPermNode();
	char ymFileName[1000], yvFileName[1000];
	sprintf(ymFileName,"YMean_ProdCondi_Iter=%d.out",iter+1);
	ofstream ym(ymFileName, ios::out);
	ym<<"##### Y Mean #####"<<endl;
	sprintf(yvFileName,"YVar_ProdCondi_Iter=%d.out",iter+1);
	ofstream yv(yvFileName, ios::out);
	yv<<"##### Y Variance #####"<<endl;

	for(int yInd1=0;yInd1<numY;yInd1++)
	{
		//Mean update
		double meann = permPtr->getYAvg(yInd1);
		for(int iMeas=0;iMeas<nMeas;iMeas++)
		{
			meann += lambda[iMeas+(yInd1*nMeas)]*prodErr[iMeas];
		}
		permPtr->setYAvg(yInd1,meann);
		ym<<yInd1<<"   "<<meann<<endl;


		//CYY update, not preserving correlation model
		for(int yInd2=yInd1;yInd2<numY;yInd2++)
		{
			double lambdaCQY = 0;
			for(int iMeas=0;iMeas<nMeas;iMeas++)
			{
				lambdaCQY += lambda[iMeas+(yInd1*nMeas)]*covQY[iMeas+(yInd2*nMeas)];
			}
			permPtr->setCYY_Cond(yInd1,yInd2,lambdaCQY);
			if(yInd1!=yInd2){permPtr->setCYY_Cond(yInd2,yInd1,lambdaCQY);}
			else{yv<<yInd1<<"   "<<permPtr->getCYY(yInd1,yInd2)<<endl;}
		}
	}
	ym.close();
	yv.close();
	permPtr->wrtPermMomnt();
}

