#include "Region.h"

Region::Region(ifstream &is, Grid *grid, Unit u ) 
 : gridPtr(grid), unit(u) {
   //cout << "Region:a:Region(). " << endl;
   x_Perm = gridPtr->getPermX();
   y_Perm = gridPtr->getPermY();
   z_Perm = gridPtr->getPermZ();

   readDeck( is );
   assignIndex();
   addCorObj();

   if(unit == FIELD) {
      convertUnit();  // from field to metric
   }

   debug = true;
   debug = false;
   if(debug) {
      wrtRegionGridPermIndex();
   }
}

Region::~Region(){
   for(int i = 0 ; i < nRegions; i++) delete corPtr[i];
   delete[] corPtr;

   delete[] yTrendType; 
   delete[] yMean; 
   delete[] yVari;
   delete[] azimuth;
   for(int i = 0; i < 3; i++) delete[] yCorr[i];
   for(int i = 0; i < 9; i++) delete[] yTrendInfo[i];
   delete[] is_perm;
   delete[] nperm;
   delete[] rInd;
   delete[] rFineInd;
   delete [] CovarianceType;
   for(int i = 0; i < 2; i++ ) {
       delete[] rx[i];
       delete[] ry[i];
       delete[] rz[i];
   }
   //interface
   for(int i = 0; i < 2; i++ ) {
       delete[] FirstRange[i];
       delete[] SecondRange[i];
   }
   delete[] InterfaceDir;
   delete[] InterfaceLoc; 
}

void Region::initialize() {
   for(int i = 0; i < 2; i++) {
       rx[i] = new int[ nRegions ]; 
       ry[i] = new int[ nRegions ]; 
       rz[i] = new int[ nRegions ];
   }
   for(int i = 0; i < 3; i++) yCorr[i]      = new double [ nRegions ];
   for(int i = 0; i < 9; i++) yTrendInfo[i] = new double [ nRegions ];
   is_perm        = new ifstream [ nRegions ];
   yMean          = new double   [ nRegions ];
   yVari          = new double   [ nRegions ];
   azimuth        = new double   [ nRegions ];
   nperm          = new int      [ nRegions ];
   CovarianceType = new int      [ nRegions ];
   yTrendType     = new int      [ nRegions ];
   corPtr         = new Covariance* [nRegions];
   rInd           = new int      [ gridPtr->getnNode() ];
   rFineInd       = new int      [ gridPtr->getNumPermNode() ]; 
   zeroingArray(rFineInd,gridPtr->getNumPermNode());//(pipatl)initialize array
}

void Region::assignIndex() {
   int ijk;
   for(int ir = 0; ir < nRegions; ir++ ) {
       for(int k = rz[0][ir]; k <= rz[1][ir]; k++ ) {
           for(int j = ry[0][ir]; j <= ry[1][ir]; j++ ) {
               for(int i = rx[0][ir]; i <= rx[1][ir]; i++ ) {
                   ijk = gridPtr->getIndex( i, j, k ); 
                   rInd[ijk] = ir;
               }
           }
       }
   }
   int k0, j0, i0, kend, jend, iend;
   for(int ir = 0; ir < nRegions; ++ir ) {
       k0 = gridPtr->toKPerm( rz[0][ir] );
       j0 = gridPtr->toJPerm( ry[0][ir] );
       i0 = gridPtr->toIPerm( rx[0][ir] );
       if(ir > 0 ) {
          if(kend < gridPtr->getPermNz() - 1) k0 = kend + 1;
          if(jend < gridPtr->getPermNy() - 1) j0 = jend + 1;
          if(iend < gridPtr->getPermNx() - 1) i0 = iend + 1;
       }
       kend = gridPtr->toKPerm( rz[1][ir] );
       jend = gridPtr->toJPerm( ry[1][ir] );
       iend = gridPtr->toIPerm( rx[1][ir] );
       for(int k = k0; k <= kend; k++ ) {
          for(int j = j0; j <= jend; j++ ) {
              for(int i = i0; i <= iend; i++ ) {
                  ijk = gridPtr->getPermIndex( i, j, k ); 
                  rFineInd[ijk] = ir;
              }
          }
       }
   }
}

//============================================================================
void Region::addCorObj() {
   for(int i = 0; i < nRegions; ++i) {
       if      ( CovarianceType[i] == 0) {
	 corPtr[i] = new ExpCovariance  (1.0, yCorr[0][i], yCorr[1][i], yCorr[2][i], azimuth[i]);
       }
       else if ( CovarianceType[i] == 1) {
	 corPtr[i] = new GaussCovariance(1.0, yCorr[0][i], yCorr[1][i], yCorr[2][i], azimuth[i]);
       }
       else {
          cerr << "no such covariance type implemented yet ";
          cerr << "please add your newe implementations" << endl;
          exit(8);
      }
   }
}

double Region::getCor(int &i1, int &j1, int &k1, int &i2, int &j2, int &k2) {
   int ijk1 = gridPtr->getPermIndex(i1, j1, k1);
   int ijk2 = gridPtr->getPermIndex(i2, j2, k2);
   if(rFineInd[ijk1] == rFineInd[ijk2] ) {
      return corPtr[rFineInd[ijk1]]->getCov(x_Perm[i1], y_Perm[j1], z_Perm[k1],
                                            x_Perm[i2], y_Perm[j2], z_Perm[k2]);
   }
   else return 0.; 
}

//============================================================================
void Region::initialize_interfaces() {
  int i;
  InterfaceDir = ( nInterfaces > 0 ) ? new int[nInterfaces] : NULL;
  InterfaceLoc = ( nInterfaces > 0 ) ? new int[nInterfaces] : NULL;
  for( i = 0; i < 2; i++ ) {
       FirstRange[i] = ( nInterfaces>0 ) ? new int[nInterfaces] : NULL;
       SecondRange[i] = ( nInterfaces>0 ) ? new int[nInterfaces] : NULL;
  }
}

// --- Read Deck ---
void Region::readDeck( ifstream &is ) {
  int i;
  Junk( is ); is >> nRegions;

  cout << endl;
  cout << "===== Region Info =====" << endl;
  cout << "Number of Regions : " << nRegions << endl << endl;
  initialize();

  int ir;
  for( ir = 0; ir < nRegions; ir++ ) {
       cout << "At Region " << ir + 1 << ":" << endl;
       Junk( is );
       is >> rx[0][ir] >> rx[1][ir]
          >> ry[0][ir] >> ry[1][ir]
          >> rz[0][ir] >> rz[1][ir];

       cout << "Range: (" << 
            rx[0][ir] << ", " << ry[0][ir] << ", " << rz[0][ir]
                           << ") to (" <<
            rx[1][ir] << ", " << ry[1][ir] << ", " << rz[1][ir]
                           << ")\n";
      
       for( i = 0; i < 2; i++ ) {
            rx[i][ir]--;
            ry[i][ir]--;
            rz[i][ir]--;
       }

       azimuth[ir] = -9999; //Default value of zero
       Junk( is );
       is >> yMean[ir] >> yVari[ir]
          >> yCorr[0][ir] >> yCorr[1][ir] >> yCorr[2][ir] >> azimuth[ir];
       if(azimuth[ir]==-9999){
	 cerr<<"Input of covariance model azimuth is required in the input deck "
	     <<"after the three principle correlation lengths"<<endl;
	 exit(8);
       }
       cout << "Permeability Mean  = " << yMean[ir] << endl;
       cout << "Permeability Vari  = " << yVari[ir] << endl;
       cout << "Permeability Corr  = " << yCorr[0][ir] <<' '
            << yCorr[1][ir] <<' '      << yCorr[2][ir] <<' '
	    << azimuth[ir] << endl;
       
       Junk(is);
       is >> yTrendType[ir]; 
            
       for( i = 0; i < 9; i++) is >> yTrendInfo[i][ir];
           
       if( yTrendType[ir] == 4 ) {
           char PermeabilityFile[1000];
           Junk(is);
           is   >> PermeabilityFile;
           cout << PermeabilityFile << endl;
           Junk(is);
           is   >> nperm[ir];
           cout << nperm[ir] << endl;
           is_perm[ir].open( PermeabilityFile, ios::in );            
       }
         
       Junk(is);
       is >> CovarianceType[ir];
       cout << "Cavriance Type = ";
       if(CovarianceType[ir] == 0) {
          cout << " Exponential ";
       } else {
          cout << " Gaussian ";
       }
       cout << " Function!" << endl;
  }
  
  // Interface
  Junk(is);
  is >> nInterfaces;
  cout <<  "nInterfaces = " << nInterfaces << endl;
  initialize_interfaces(); // initialize interface scope data
  int iface;
  for(iface = 0; iface < nInterfaces; iface++ ) {
      Junk(is);
      is >> InterfaceDir[iface] >> InterfaceLoc[iface]
         >> FirstRange[0][iface] >> FirstRange[1][iface]
         >> SecondRange[0][iface] >> SecondRange[1][iface]; 
                   
      // --- back to zero based index ---
      InterfaceLoc[iface]--;
                                  
      for(i = 0; i < 2; i++ ) {
          FirstRange[i][iface]--;
          SecondRange[i][iface]--;
      }
   }
}

void Region::wrtRegionGridPermIndex() {
  char outfile[20];
  ofstream os("RegnInfo.out", ios::out);
  os << setiosflags(ios::fixed | ios::showpoint)
     << setprecision(4);
 
  int ijk;
  for(int k = 0; k < gridPtr->getPermNz(); ++k) {
     for(int j = 0; j < gridPtr->getPermNy(); ++j) {
         for(int i = 0; i < gridPtr->getPermNx(); ++i) {
             ijk = gridPtr->getPermIndex(i, j, k);
             os << x_Perm[i] << ' ' 
                << y_Perm[j] << ' '
                << z_Perm[k] << ' '
                << rFineInd[ijk]<<' ' 
                << endl;
         } 
         os << endl;
     }
  }
  os.close();
}

//===== convertUnit =====
void Region::convertUnit() {
   for(int j = 0; j < 3; j++) {
      for(int i = 0; i < nRegions; i++) yCorr[j][i] *= ft_m;
   }
}

//(pipatl)To initialize new array(primarily for rFineInd)
void Region::zeroingArray(int *array, int size)
{
	for(int ii=0;ii<size;ii++)
	{
		array[ii] = 0;
	}
}
