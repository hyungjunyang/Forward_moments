#ifndef _QUADRATURE_H_
#define _QUADRATURE_H_

#include <iostream>
#include <math.h>

using namespace std;


void GaussianQuadrature(int qPoint, double* &rloc, double* &weight)
{
	rloc = new double[qPoint];
	weight = new double[qPoint];
	switch(qPoint)
	{
	        case 1:
		        rloc[0] = 0;
			weight[0] = 2;
			break;
		case 2:
			rloc[0] = -sqrt(1.0/3.0);
			rloc[1] = sqrt(1.0/3.0);
			weight[0] = weight[1] = 1;
			break;	
		case 3:
			rloc[0] = -sqrt(3.0/5.0);
			rloc[1] = 0;
			rloc[2] = sqrt(3.0/5.0);
			weight[0] = weight[2] = 5.0/9.0;
			weight[1] = 8.0/9.0;
			break;
		case 5:
			rloc[0] = -sqrt(5.0+(2.0*sqrt(10.0/7.0)))/3.0;
			rloc[1] = -sqrt(5.0-(2.0*sqrt(10.0/7.0)))/3.0;
			rloc[2] = 0;
			rloc[3] = sqrt(5.0-(2.0*sqrt(10.0/7.0)))/3.0;
			rloc[4] = sqrt(5.0+(2.0*sqrt(10.0/7.0)))/3.0;
			weight[0] = weight[4] = (322.0-(13.0*sqrt(70.0)))/900.0;
			weight[1] = weight[3] = (322.0+(13.0*sqrt(70.0)))/900.0;
			weight[2] = 128.0/225.0;
			break;
	        case 7:
		        rloc[0] = -0.94910791;
			rloc[1] = -0.74153119;
			rloc[2] = -0.40584515;
			rloc[3] = 0;
			rloc[4] = 0.40584515;
			rloc[5] = 0.74153119;
			rloc[6] = 0.94910791;
			weight[0] = 0.12948497;
			weight[1] = 0.27970539;
			weight[2] = 0.38183005;
			weight[3] = 0.41795918;
			weight[4] = 0.38183005;
			weight[5] = 0.27970539;
			weight[6] = 0.12948497;
			break;
	        case 10:
		        rloc[0] = -0.97390653;
		        rloc[1] = -0.86506337;
		        rloc[2] = -0.67940957;
		        rloc[3] = -0.43339539;
		        rloc[4] = -0.14887434;
		        rloc[5] = 0.14887434;
		        rloc[6] = 0.43339539;
		        rloc[7] = 0.67940957;
		        rloc[8] = 0.86506337;
		        rloc[9] = 0.97390653;
			weight[0] = 0.06667134;
			weight[1] = 0.14945135;
			weight[2] = 0.21908636;
			weight[3] = 0.26926672;
			weight[4] = 0.29552422;
			weight[5] = 0.29552422;
			weight[6] = 0.26926672;
			weight[7] = 0.21908636;
			weight[8] = 0.14945135;
			weight[9] = 0.06667134;
			break;
		default:
			cerr<<qPoint<<"-point quadrature rule hasn't been implemented yet!!!"<<endl<<endl;
	}
}

#endif

