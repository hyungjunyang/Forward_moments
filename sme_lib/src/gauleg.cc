//////////////////////////////////////////////////
//                                              //
//                 Liyong Li                    //
//      Reservoir Simulation Research Team      //
//        Chevron Petroleum Technology Co.      //
//          1300 Beach Blvd., Rm. 3166          //
//            La Habra, CA 90631-6374           //
//             Phone:(562) 694-7366             //
//              Fax:(562) 694-7565              //
//            Email liyl@chevron.com            //
//                                              //
//                 Version 1.                   //
//                Apr. 07, 1999                 //
//////////////////////////////////////////////////

#include <math.h>
#include <iostream>
#define EPS 3.0e-11

void gauleg(double x1, double x2, double x[], double w[], int n)
{
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;
       
//      adding int k for shift array band from [1,n] to [0,n-1]

        int k;
 
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);

//		x[i]=xm-xl*z;
//		x[n+1-i]=xm+xl*z;
//		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
//		w[n+1-i]=w[i];
//                cout<<i<<' '<<x[i]<<' '<<w[i]<<endl;

                k=i-1;
                x[k]=xm-xl*z;
                x[n-1-k]=xm+xl*z;
                w[k]=2.0*xl/((1.0-z*z)*pp*pp);
                w[n-1-k]=w[k];

//                cout<<k<<' '<<x[k]<<' '<<w[k]<<endl;

	}
}
#undef EPS
