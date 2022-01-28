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

//          --       Domain.h          --
//          --    class definition    --

#ifndef DOMAINH
#define DOMAINH
#include <iostream>
#include <fstream>

#include "Common.h"
#include "Control.h"
#include "Grid.h"
#include "Region.h"
#include "Fluid.h"

#include "Block.h"
#include "WellArc.h"
#include "Well.h"

class Domain{

public:
   Domain(ifstream &is, char *Dir);
   Domain(Grid*, int, Well* Well);
   virtual ~Domain();	// destructor
   
   void addBlocks();
   //void buildBlocks();
   //void buildBlocks(int, int, double, double, double, double, double, double, double, double);
   //void buildWellArcs(int, int, bool, bool, double, double, double, double, double);
   void addWellArcs();

   // read()
   void readMCVmoments();

   // add()
   void addWellObj(int, int, bool);
   void addVeloCova(double*, double*, double*);
   
   // delete() // (Pipatl) This is to be used before calling Flow::writeToDomain() repeatedly
   void deleteVeloCovariances();

   // set()
   void setVelocity(double*, double*);

   void setBoundBlk(); 
   // calc()
   void calcDerivative();

   // get()
   Grid*    getGrid()   {return gridPtr;}
   Fluid*   getFluid()  {return fluidPtr;}
   Region*  getRegn()   {return regnPtr;}
   int      getDtNum()  {return nDt;}
   
   Control* getContr()        {return contPtr;}
   double getSwc()            {return contPtr->getSwc();}
   double getSor()            {return contPtr->getSor();}
   double getViscosR()        {return contPtr->getViscosR();}
   double getProductionTime() {return contPtr->getProductionTime();}
   int    getProductionTimeStep() {return contPtr->getProductionTimeStep();}

   int      getWellNum(){return nWells;}
   Well*    getWellArr(){return wellArr;}  
   
   double   getVolume() {return volume;}
   WellArc* getWellArc(int i) {return wellArcArray[ i ];}
   Block * getBlock(int i, int j, int k) {return blocks[ getIndex(i,j,k) ];}
   Block * getBlock(int i, int j) {return blocks[ getIndex(i,j) ];}
   Block * getBlock(int i) {return blocks[ i ];}
   //Block * getProdBlock()  {return prodBlk; }

   int getNx() {return nx_;}
   int getNy() {return ny_;}
   int getNz() {return nz_;}
   Unit getUnit() {return unit;}
   double getBdx(int i) {return gridPtr->getBdx(i);}
   double getBdy(int i) {return gridPtr->getBdy(i);}
   double getBdz(int i) {return gridPtr->getBdz(i);}
   int getIndex(int i, int j) {return i + j * nx_;}
   int getIndex(int i, int j, int k) {return i+(j*nx_)+(k*nx_*ny_);}
   int getNumBlocks() {return node_;}
   int getNumWells()  {return nWells;}
   double* getCv1v1()    {return Cv1v1;}
   double* getCv2v2()    {return Cv2v2;}
   double* getCv1v2()    {return Cv1v2;}
   double* getDxCv1v1()  {return dxCv1v1;}
   double* getDyCv1v1()  {return dyCv1v1;}
   double* getDxCv2v2()  {return dxCv2v2;}
   double* getDyCv2v2()  {return dyCv2v2;}
   double* getDx1Cv1v2() {return dx1Cv1v2;}
   double* getDy1Cv1v2() {return dy1Cv1v2;}
   double* getDx2Cv1v2() {return dx2Cv1v2;}
   double* getDy2Cv1v2() {return dy2Cv1v2;}
   double  getPoro(int i, int j, int k) {
	           return fluidPtr->getPoro( gridPtr->getIndex(i,j,k)
			             ); 
   }
   // misc
   Block* point2block(Point& b);
   void changeBlockVeloSign();
   void changeBlockProdInjc();
   void  calcVolume();
   bool IsInDomain(int i,int j,int k){return (i>=0)&&(i<nx_)&&(j>=0)&&(j<ny_)&&(k>=0)&&(k<nz_);}

   //print()
   void PrintDomainInfo(ostream & os);
   void printWellArc();
private:
   bool debug;
   
   Unit unit;
   void readUnit (ifstream &is );
   
   Grid *gridPtr;
   void addGridObj(ifstream &is, char * Dir);

   Region *regnPtr;
   void addRegnObj(ifstream &is, Unit unit);

   Fluid *fluidPtr;
   void addFluidObj(ifstream &is, Unit unit );
  
   int nWells;
   Well  *wellArr;
   WellArc** wellArcArray;
   void  addWellArr(ifstream &is );

   int nDt;
   Control *contPtr;
   void  addContObj(ifstream &is );
   

   int nx_, ny_, nz_, node_;
   double volume;
   Block** blocks;
   Block* wellBlk;
   Block* boundBlk; 
   double *Cv1v1, *Cv2v2, *Cv1v2;
   double *dxCv1v1,  *dyCv1v1,  *dxCv2v2,  *dyCv2v2; 
   double *dx1Cv1v2, *dx2Cv1v2, *dy1Cv1v2, *dy2Cv1v2;

};

#endif 

