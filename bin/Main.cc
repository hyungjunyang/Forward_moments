/*
 * File: Main.cxx
 * ----------------------------------------------------------------------
 * The main program for whole program (stochastic formulation). This Main 
 * program includes two separate programs, one for flow calculation, the 
 * other for transport calculation. They will run separately to 
 * conserve memory.  
 * 
 * The connection between these two programs: flow calculation will output 
 * several files containing grid geometry and velocity moments, transport 
 * calcultion will read these files, and start from there.
 *
 * If you want, you can separate this file into two small files, one for 
 * flow, the other for transport. A quite simple process.
 * 
 * ----------------------------------------------------------------------
 * Ways to run this program (suppose its name is 'sim')
 *    sim 0:  flow calculation only 
 *    sim 1:  transport calculation only
 *    sim 2:  both calculation
 *    sim  :  remind you to choose the calculation
 * Note: if the argv (or choice) is not 0,1or2, you have chance to input again
 *
 * Each calcultion requires a input file, program will remind you to input
 * one, if not entered (just return), program will use the default one:
 *    flow calculation:        "flow.in"  is the default input
 *    transport calculation:   "trans.in" is the default input
 *
 * ----------------------------------------------------------------------
 * This program run in two unit,  Field or Metric, which will be used 
 * for input and output purpose only.
 * The internal calculation unit is always metric, make everything consistent.
 * Field Unit: [length]=ft, [time]=day, [K=k/mu]=md/cp, [p] = psia, 
 *             [cf,cr]=1/psia, [v]=ft/day, [Q]=STB/day, [den]=lbm/SCF
 *
 * ----------------------------------------------------------------------
 *
 */

#include "Flow.h"
#include "ParticleTrack.h"
#include "Saturation.h"
#include "MC2SME.h"

int main(int argc, char *argv[]) {
  int choice;
  char dummy[1000], inputFile1[1000], inputFile2[1000], directory[1000];
  ifstream is_flow, is_trans;

  // --- input the choice for calculation (if necessary) -----
  while(true) {
    if(argc > 1) { 
       choice = atoi(argv[1]);
    } else {
       cerr << endl;
       cerr << "********Calculation choice:*****************" << endl;
       cerr << "**   0: flow calculation only             **" << endl;
       cerr << "**   1: flow and travel time calculations **" << endl;
       cerr << "********************************************" << endl << endl;
       cerr << "Input Calculation choice: ";
       cin >> choice;
       cin.getline(inputFile1, 1000);  // eat the "\n"
    }
    if(choice >=0 && choice <= 2) break;
    argc = -1;
    cerr << "\n" <<choice << " is a WRONG choice, input again!!! " << endl;
  }
    
  cerr << "Enter working directory: ";
  cin.getline(dummy, 1000);
  if(dummy[0] != '\0') sprintf( directory, "%s/", dummy );
  else sprintf( directory, "%s", dummy );

  // --- read the input file names -------- 
  cerr << "Flow part Input file name (flow.in): ";
  cin.getline(dummy, 1000);
  if(dummy[0] != '\0')
     sprintf(  inputFile1, "%s%s", directory, dummy );
  else
     sprintf(  inputFile1, "%sflow.in", directory );

  is_flow.open(inputFile1,ios::in);
  if( is_flow.bad() ) {
      cerr << " Can not open " <<  inputFile1 << endl;
      exit(8);
  }
  bool debug = false;
  if( debug ) {
     cout << inputFile1 << endl;
  }
  
  // ===== Module One =====
  Domain *domainPtr = new Domain(is_flow, directory);
  /*
  cout << domainPtr->getSwc() << ' '
       << domainPtr->getSor() << ' '
       << domainPtr->getViscosR() << ' '
       << domainPtr->getProductionTime() <<' '
       << domainPtr->getProductionTimeStep()<<' '
       << endl;
  */
  // ===== Module Two =====
  Flow *flowPtr = new Flow(directory, domainPtr, debug);

  flowPtr->solve(directory, is_flow);
  flowPtr->writeToFile(directory);

  if( choice == 0) {
	  delete flowPtr;
	  delete domainPtr;
      cout<<"Calculations of Flow Moments Are Completed!!";
      cout<<endl<<endl;
      return 0;
  }
  flowPtr->writeToDomain( domainPtr );

  // ===== Module Two =====
  double startT = clock();
  cout<<endl;
  cout<<"===== Starting Calculations of Travel Time Moments =====" <<endl;
          
  //domainPtr->readMCVmoments();

  ofstream fileOu("block.out", ios::out);
  //domainPtr->changeBlockVeloSign();
  //domainPtr->changeBlockProdInjc();
  domainPtr->PrintDomainInfo(fileOu); 
  fileOu.close();
  
  // -- start streamline tracking
  //debug = true;
  ParticleTrack* ptclTrackPtr = new ParticleTrack(domainPtr, debug, true);
  ptclTrackPtr->goTrack();
  
  cout<<endl;
  cout<<"===== The Time for Calculations of Travel Time Moments Are Completed are: ";
  double endT = clock();
  cout << (endT-startT)/CLOCKS_PER_SEC 
       << " sec =====" << endl;

  
  delete flowPtr;
  delete ptclTrackPtr;
  delete domainPtr;


  cout << endl;
  cout << "The program is completed normally!" << endl;
  cout << endl;
  
  return 0;
}

