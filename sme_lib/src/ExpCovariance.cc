#include "ExpCovariance.h"

ExpCovariance::ExpCovariance(double var, double cor1, double cor2, double cor3, double azimuth)
: Covariance (var, cor1, cor2, cor3, azimuth)
{
  cout << "ExpCovariance::ExpCovariance() \n ";
}

ExpCovariance::~ExpCovariance(){
}

double ExpCovariance::getCov(double x1, double y1, double z1, 
                             double x2, double y2, double z2){
  double hxr,hyr,hzr; //rotated coordinate of the distances
  hxr = ((x1-x2)*cos(azimuth*PI/180))+((y1-y2)*sin(azimuth*PI/180));
  hyr = ((x1-x2)*-1*sin(azimuth*PI/180))+((y1-y2)*cos(azimuth*PI/180));
  hzr = z1-z2;

  double sum = 0.;
//cout << "ExpCovariance::getCov()"<<endl;
  sum  = pow( hxr / cor_length1, 2. );
  sum += pow( hyr / cor_length2, 2. );
  sum += pow( hzr / cor_length3, 2. );
  sum =  variance * exp( -sqrt( sum ) );
  return sum;
}
