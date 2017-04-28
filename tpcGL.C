/*
  root -q -b lGarfield.C tpcGL.C+
*/

#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TString.h"
#include "ViewField.hh"
#include "ViewDrift.hh"
#include "ViewCell.hh"
#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentAnalyticField.hh"
#include "Sensor.hh"
#include "DriftLineRKF.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"

#include "FundamentalConstants.hh"
#include "Random.hh"
#include "Plotting.hh"
//#define __DriftVelocity__
using namespace Garfield;
TFile *fOut = 0;
//________________________________________________________________________________
void verticalWall(ComponentAnalyticField* comp, double x, double yMin, double yMax, double voltage, const char * name)
	{
//	FILE * wallCoords;
//	wallCoords = fopen("verticalWAll.txt", "a");
	double step = 0.011;
	double dPseudoWire = 0.01;

	int numberOfWires = (yMax - yMin) / step;
	for (double y = yMin; y < yMax; y += step)
		{
		comp->AddWire(x, y, dPseudoWire, voltage, name);
//		fprintf (wallCoords, "%s - %f - %f - %f\n", name, x, y, dPseudoWire);
		}
//	fclose(wallCoords);
	}
void horizontalWall(ComponentAnalyticField* comp, double xMin, double xMax, double y, double voltage, const char * name)
	{
	double step = 0.011;
	double dPseudoWire = 0.01;

	int numberOfWires = (xMax - xMin) / step;
	for (double x = xMin; x < xMax; x += step)
		comp->AddWire(x, y, dPseudoWire, voltage, name);
	}

void TPCAnodeWires(ComponentAnalyticField* comp, const Char_t *setup, 
		   double xmin, double xmax, 
		   Int_t na,            /* total no. of anode wires */   
		   Int_t nf,         	/* no. of field wires       */   
		   double xW,      	/* X position of zero wire */    
		   double yW,      	/* Y position of zero wire */    
		   double aWpitch, 	/* anode wire pitch */	      
		   double da,      	/* anode sense wire diamiter */  
		   double df,      	/* field wire diameter */	      
		   double V,	     	/* Voltage */                    
		   Bool_t zFLW = kFALSE /* kTRUE if Zero field last wire */
		   )
	{
	cout << "TPCAnodeWires: " << setup << " total no. of anode wires = " << na + 1
		<< " no. of field wires from each side " << nf 
		<< " X position of zero wire " << xW 
		<< " Y position of zero wire " << yW 
		<< " anode wire pitch " << aWpitch
		<< " anode sense wire diamiter " << da 
		<< " field wire diameter " << df 
		<< " Voltage " << V
		<< " kTRUE if Zero field last wire " << zFLW << endl;
	const Char_t *s[3] = {"s","S","Z"};
	double x;
	Int_t c = 0;

//	FILE * anodeCoords;
//	anodeCoords = fopen("anodeCoords.txt", "a");

	for (Int_t i = 0; i <= na; ++i)
		{
		x = xW + aWpitch*i;
		if (x < xmin || x > xmax) continue;
		double d = da;
		double volt = V;
		c = 0;
		if (i == na) continue; // removing last fat wire
		if (i == na - 1) continue; // removing second to last fat wire
		if (i < nf || i > na - nf)
			{
			d = df;
			c = 1;
			}
		if (zFLW && (i == 0 || i == na))
			{
			volt = 0;
			c = 2;
			}
		comp->AddWire(x, yW, d, volt, s[c]);
//		fprintf(anodeCoords, "%f\t%f\t%f\n", x, yW, d);
		continue;
		}
//	fclose(anodeCoords);
	}
//________________________________________________________________________________
void tpcGL(Int_t nEvents = 0, const Char_t *OutName = "GL.root", TString geoName = "PWD")
	{
	TString Geometry;
	if (geoName == "PWD")
		{
		TString PWD(gSystem->WorkingDirectory());
		Geometry = (TString) gSystem->BaseName(PWD);
		cout << "Directory : " << Geometry.Data() << endl;
		}
	else Geometry = geoName;
	TString Region("IO");
	if (Geometry.Contains("II"))      Region = "II";
	else if (Geometry.Contains("OO")) Region = "OO";  
	cout << "Region\t" << Region.Data() << endl;

	if (! (	Geometry.BeginsWith("TPC") ||
		Geometry.BeginsWith("myGeo") ||
		Geometry.BeginsWith("iTPC_3x125mkm") ||
		Geometry.BeginsWith("iTPC_3x250mkm") ||
		Geometry.BeginsWith("iTPC_3x500mkm") ||
		Geometry.BeginsWith("iTPC_3x1mm")    ||    
		Geometry.BeginsWith("iTPC_3x125mkm_JT_091515") ||
		Geometry.BeginsWith("iTPC_3x125mkm_JT_091715") ||
		Geometry.BeginsWith("iTPC_3x125mkm_JT_092215_Proposal") 
		))
		{
		cout << "Geometry has not been recognized" << endl;
		return;
		}

	cout << "Geometry\t" << Geometry.Data() << endl;
//	STAR coordinate system (xS,yS,zS) => Garfield (yS,zG,xG,yG); Garfield(xG,yG,zG) = > Star(y,z,x)
//	plottingEngine.SetDefaultStyle();
//	Setup the gas.
	const double BarPressure         = 1010.8; // [mbar], TPC-PTB, barometricPressure
	const double inputTPCGasPressure = 1.93;   // [mbar], TPC-PT8, difference between barometer pressure and pressure in TPC
//	const double pressure = (1011. / 1000.) * 750.; // 1 bar = 750 torr 
	const double pressure = ((BarPressure + inputTPCGasPressure)/ 1000.) * 750.; // [torr], 1 bar = 750 torr 
	const double temperature = 297.839; //273.15 + 24.7; // inputGasTemperature (degree K)
	double BField = 0.5; // Tesla
	double Angle  = 0.0; // rad
	TString gasFile("P10.");
	Int_t b = 10*BField;
	Int_t t = temperature;
	Int_t p = pressure;
	gasFile += Form("B%ikGT%iP%i",b,t,p);
//	Switch between Inner and Outer sector
//	Voltage settings
	double vAnodeI  =   1100.;
	if (Geometry.Contains("1170")) vAnodeI = 1170.;
	const double vAnodeO  =   1390.;
//	Create a gas medium.  
	MediumMagboltz* gas = new MediumMagboltz();
#if 0
	const double density = LoschmidtNumber * 
		(pressure / AtmosphericPressure) *
		(273.15 / temperature);
	const double pt = (temperature / 273.15) * (760. / pressure);
#endif
//	Set the temperature [K] and pressure [Torr].
	gas->SetTemperature(temperature);
	gas->SetPressure(pressure);
	gas->SetMaxElectronEnergy(300.);
//	Specify the gas mixture.
	double Fracs[2] = {0.9, 0.1};
	gas->SetComposition("ar", 100*Fracs[0], "ch4", 100*Fracs[1]);
	gas->SetFieldGrid(100.,100e3, 20, true, BField, BField, 1, Angle, Angle, 1);
	gasFile += ".gas";
	const double rPenning = 0.57;
	const double lambdaPenning = 0.e-4;
	gas->EnablePenningTransfer(rPenning, lambdaPenning, "ch4");
	if (gSystem->AccessPathName(gasFile))
		{
//		Use Magboltz to generate the gas table and grid.
//		Specfify te number of collission [1e7] over which the electron is traced in Magboltz
		const Int_t ncoll = 10; //[1e7]
		const Bool_t verbose = true;
		gas->GenerateGasTable(ncoll, verbose);
//		Save the gas table for later use
		gas->WriteGasFile(gasFile.Data());
		return;
		}
	else
		{
		gas->LoadGasFile(gasFile.Data());
		}
//	gas->PrintGas();
	const Char_t *path = gSystem->ExpandPathName("$GARFIELD_HOME/Data/IonMobility_Ar+_Ar.txt");
	gas->LoadIonMobility(path);
	delete path;
	gas->Initialise();
//	Geometry/tpc/tpcWirePlanes
	double  anodeWireRadius                      = 0.0020/2;// 20 um diameter 
	double  frischGridWireRadius                 = 0.0075/2;// 
	double  gatingGridWireRadius                 = 0.0075/2;// 
	double  anodeWirePitch                       = 0.4;// 
	double  frischGridWirePitch                  = 0.1;// 
	double  gatingGridWirePitch                  = 0.1;// 
	double  innerSectorAnodeWirePadSep           = 0.2;// AnodeWire-to-PadPlane distance  
	double  innerSectorFrischGridPadSep          = 0.4;// FrischGrid-to-PadPlane distance  
	double  innerSectorGatingGridPadSep            = 1;// GatingGrid-to-PadPlane distance  
	double  outerSectorAnodeWirePadSep           = 0.4;// AnodeWire-to-PadPlane distance  
	double  outerSectorFrischGridPadSep          = 0.8;// FrischGrid-to-PadPlane distance  
	double  outerSectorGatingGridPadSep          = 1.4;// GatingGrid-to-PadPlane distance  
	Int_t     numInnerSectorAnodeWires             = 170;// 
	Int_t     numInnerSectorFrischGridWires        = 681;// 
	Int_t     numInnerSectorGatingGridWires        = 681;// 
	double  firstInnerSectorAnodeWire            = 53.2;// 
	double  firstInnerSectorFrischGridWire       = 53;// 
	double  firstInnerSectorGatingGridWire       = 53;// 
	double  lastInnerSectorAnodeWire             = 120.8;// 
	Int_t     numOuterSectorAnodeWires             = 172;// 
	Int_t     numOuterSectorFrischGridWires        = 689;// 
	Int_t     numOuterSectorGatingGridWires        = 689;// 
	double  firstOuterSectorAnodeWire            = 122.8;// 
	double  firstOuterSectorFrischGridWire       = 122.59;// 
	double  firstOuterSectorGatingGridWire       = 122.59;// 
	double  lastOuterSectorAnodeWire             = 191.19;// 
	double  fieldWireRadius                      = 0.0125/2;

//	TPC Sector 12  => Garfield : (x pad,y raw, z drift) => ( y, z, x);
//	Drawing 24A055B
	const double xMinInnerStrongBack = 51.905;
	const double xMaxInnerStrongBack = xMinInnerStrongBack + 69.527;
	const double xMinInnerGGWire     = 53.00;
	const double xMaxInnerGGWire     = xMinInnerGGWire + 68.0;
	const double xMinOuterStrongBack = 121.732;
	const double xMaxOuterStrongBack = xMinOuterStrongBack + 71.514;
	const double xMinOuterGGWire     = 122.595;
	const double xMaxOuterGGWire     = xMinOuterGGWire + 68.8;

	const double yGG      =      0.; // Gating Grid
	const double vGG      =   -115.; // Gating Grid
	const double dvGG     =     75.; // alternative delta Voltage for closed GG
	const double vCathode =  -27950;
	const double yCathode =-208.707; // move into system where yGG = 0
	const double eDrift   = (vCathode - vGG)/(yCathode - yGG); 
	const double yC       = -1.0;     // pseudo cathode
	const double vC       = vGG + eDrift*(yC - yGG);
	const double gapI     = 0.2; // Anode pad plane gap Inner
	const double gapO     = 0.4; // Anode pad plane gap Outer
	const double aWpitch  = 0.4; // Anode wire spacing;
	const double cWpitch  = 0.1; // Cathode wire spacing
	const double ggWpitch = 0.1; // 
//	y coordinates of the wires
	const double yFG      = 0.6;
	const double vPad     = 0;

//	Periodicities
  
//	Wire diameters
	const double dSens = 2*anodeWireRadius;
	const double dFG   = 2*frischGridWireRadius;
	const double dGG   = 2*gatingGridWireRadius;
	const double dFW   = 2*fieldWireRadius; // Field Wire
	const Int_t nBinsAngular = 360;
	const Int_t nBinsRadial = 100;
	const Int_t nBinsGain = 100;
	double ymin = yC; 
	double ymax = yGG + yFG + 2*gapO;
	double dZ   = 100;
	double zmin = -dZ, zmax = dZ;
	double xmin, xmax;
	if (Region == "IO")
		{
		xmin = 0.0100*(TMath::Nint((xMaxInnerGGWire - 3.0)/0.0100) - 1); // Inner to Outer sector transition region 
		xmax = 0.0100*(TMath::Nint((xMinOuterGGWire + 3.0)/0.0100) + 1); //
		}
	else if (Region == "II")
		{
		xmin = 0.0100*(TMath::Nint((xMinInnerGGWire - 3.0)/0.0100) - 1); // Inner 
		xmax = 0.0100*(TMath::Nint((xMinInnerGGWire + 3.0)/0.0100) + 1); //
		}
	else if (Region == "OO")
		{
		xmin = 0.0100*(TMath::Nint((xMaxOuterGGWire - 3.0)/0.0100) - 1); // Outer
		xmax = 0.0100*(TMath::Nint((xMaxOuterGGWire + 3.0)/0.0100) + 1); //
		}
/*
y
^ 
-> x
       Inner                                Outer
       
                                        | | ____________________________             yPad[1]
                                        | | xSBmin[1]                   xSBmax[1]
					| |                 
					| |
					| |
					| |              gapO
					| |
					| |
					| |
                             dXWall	| |
yPad[0] ________________________| |	| | ..........................    AnodeW ySens[1], 
    xSBmin[0]          xSBmax[0]| |     | |
                                | |     | |
       gapI                     | |     | |
                                | |     | |
      ........................  | |     | |                               AnodeW ySens[0]
                                | |     | |
                                | |     | |
                                | |     | |
                                | |     | |
      ************************  |-|     |-|-- **************************    yFG
                                         




      ........................               ..........................    yGG = 0
      xGGmin[0]                              xGGmin[1]

      ------------------------------------------------------------- Cathode  
 







*/
//	Build the geometry, in this case just a box.
	GeometrySimple* geo = new GeometrySimple();
//	geo->EnableDebugging();
	SolidBox* box = new SolidBox(0.5*(xmax + xmin), 0.5*(ymax + ymin), 0.5*(zmax + zmin), 
			       0.5*(xmax - xmin), 0.5*(ymax - ymin), 0.5*(zmax - zmin));
//	Add the solid to the geometry, together with the medium inside.
	geo->AddSolid(box, gas);

//	Setup the electric field.
	ComponentAnalyticField* comp = new ComponentAnalyticField(); 
	comp->SetMagneticField(0., BField, 0.);
//	comp->SetPeriodicityX(nRep * perSens);
 
	double dXWall = 0.18; // Jim Thomas wall 
	double VWall  = -115;
	double VWall2 = 0;
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_091515_200V")) VWall2 = -200;
	if (Geometry.BeginsWith("myGeo")) VWall2 = -200;
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_091715_230V")) VWall  = -230;
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_091715_460V")) VWall  = -460;
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_091715_690V")) VWall  = -690;
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_092215_Proposal")) VWall = 0;
	cout << "VWall = " << VWall << "\tVWall2 = " << VWall2 << endl;
	double xGGmin[2] = {xMinInnerGGWire, xMinOuterGGWire};
	double xGGmax[2] = {xMaxInnerGGWire, xMaxOuterGGWire};
	double vAnode[2] = {vAnodeI, vAnodeO};
	double xSBmin[2] = {xMinInnerStrongBack, xMinOuterStrongBack};
	double xSBmax[2] = {xMaxInnerStrongBack, xMaxOuterStrongBack};
	if (Geometry.BeginsWith("iTPC_3x125mkm_JT_092215_Proposal"))
		{
		cout << "Correct xSBmax[0] = " << xSBmax[0] << " and xSBmin[1] = " << xSBmin[1] << " by dXWall = " << dXWall << endl;
		xSBmax[0] -= dXWall;
		xSBmin[1] += dXWall; 
		cout << " new xSBmax[0] = " << xSBmax[0] << " and xSBmin[1] = " << xSBmin[1] << endl;
		}
	else if (Geometry.BeginsWith("iTPC_3x125mkm_JT_09"))
		{
		cout << "Correct xSBmax[0] = " << xSBmax[0] << " by dXWall = " << dXWall;
		xSBmax[0] -= dXWall;
		cout << " new xSBmax[0] = " << xSBmax[0] << endl;
		}
	Int_t        ng[2] = {(Int_t) ((xMaxInnerGGWire - xMinInnerGGWire)/ggWpitch), (Int_t) ((xMaxOuterGGWire - xMinOuterGGWire)/ggWpitch)};
	Int_t        na[2];
	double    gap[2] = {gapI, gapO};
	double  ySens[2] = {yFG + gap[0], yFG + gap[1]};
	double   yPad[2] = {ySens[0] + gap[0], ySens[1] + gap[1]};
	double x;

//	FILE * coordsFile;
//	coordsFile = fopen("coords.txt", "w");

	for (Int_t io = 0; io < 2; io++)
		{
//		Inner/Outer Loop
//		Add the Gating Grid wires.
		for (Int_t i = 0; i <= ng[io]; ++i)
			{
			x = xGGmin[io] + ggWpitch*i;
//			if (!io) fprintf(coordsFile, "%f\n", x);
			if (x < xmin || x > xmax) continue;
//			double v = vGG + dvGG*(1 - 2*((i + io)%2)); // GG Closed
			double v = vGG + dvGG*(1 - 2*((i + 1 + io)%2)); // GG Closed // + 1 to shift the alternating potential by one wire
//			double v = vGG; // GG Open
			if (io == 0 && i == ng[io]) {cout << "The last  GG wire Inner Sector potential " << v << endl;}
			if (io == 1 && i ==      0) {cout << "The first GG wire Outer Sector potential " << v << endl;}
			comp->AddWire(x, yGG, dGG, v, "g");
//			Add the cathode (frischGrid) wires
			comp->AddWire(x, yFG, dFG, 0., "c");
			}
//		Add the anode wires.
		na[io] = ng[io]*(ggWpitch/aWpitch) - 1;
		if (io)
			{
//			Outer
			Int_t nf = 1;
			TPCAnodeWires(comp, "Outer", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, dFW, vAnode[io]);
			}
		else
			{
//			Inner
			Bool_t zFLW = kFALSE;
			if (Geometry.Contains("ZFLW")) zFLW = kTRUE;
			if (Geometry.BeginsWith("TPC"))
				{
				Int_t nf = 1;
				TPCAnodeWires(comp, "Inner", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, dFW, vAnode[io], zFLW);
				}
			else if (Geometry.BeginsWith("iTPC_3x125mkm"))
				{
				Int_t nf = 3;
				TPCAnodeWires(comp, "Inner", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, dFW, vAnode[io], zFLW);
				}
			else if (Geometry.BeginsWith("iTPC_3x250mkm"))
				{
				Int_t nf = 3;
				TPCAnodeWires(comp, "Inner", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, 2*dFW, vAnode[io], zFLW);
				}
			else if (Geometry.BeginsWith("iTPC_3x500mkm"))
				{
				Int_t nf = 3;
				TPCAnodeWires(comp, "Inner", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, 4*dFW, vAnode[io], zFLW);
				}
			else if (Geometry.BeginsWith("myGeo"))
				{
				Int_t nf = 3; // Number of field wires - Fat Anode Wires
				TPCAnodeWires(comp, "Inner", xmin, xmax, na[io], nf, xGGmin[io] + ggWpitch*2, ySens[io], aWpitch, dSens, dFW, vAnode[io], zFLW);
				}

			}
//		Replace pad plane by wires for inner sector
		if (!io)
			{
			double dPseudo = 0.0100; // diameter of pseudo pad wire 
			Int_t np = (xSBmax[io] - xSBmin[io])/(2*dPseudo);
			for (Int_t i = 0; i <= np; i++)
				{
				x = xSBmin[io] + dPseudo*(2*i + 0.5);
				if (x < xmin || x > xmax) continue;
				comp->AddWire(x, yPad[io], dPseudo, vPad, "P");
//				fprintf(coordsFile, "fake pad : %f - %f - %f\n", x, yPad[io], dPseudo);
				}
			}
		else
			{
			comp->AddPlaneY(yPad[io], vPad, "p");
			}
		} // Inner / Outer Loop
//	fclose(coordsFile);

	double Xmin = xGGmin[0] + ggWpitch*2 + ggWpitch*ng[0] + 0.2; // 2 mm safety margin
	double Xmax = xGGmin[1] - ggWpitch*2 - 0.2;
	double Xmean= 0.5*(Xmin + Xmax);
	double xstep = 0.0100;

//	bring wall of wires down to the level of the outer padplane
	double padStep = 0.011;
	verticalWall(comp, xSBmax[0] + padStep, yPad[0], yPad[1] - padStep / 2, vPad, "innerPadToOuterPadRadius");
	comp->AddReadout("innerPadToOuterPadRadius");

//	======================= MY GEOMETRY ================================================
	if (Geometry.BeginsWith("myGeo_V01"))
		{
		double step = 0.011;

//		Inner Wall
		double myV = -40; 

//		position at myV
		double myY = yFG - (myV / vGG) * yFG;

//		y positions of the walls
		double topY = myY + dXWall / 2 - step / 2;
		double lowY = myY - dXWall / 2 + step / 2;

		verticalWall(comp, xSBmax[0] - dXWall + step / 2,	topY, yPad[0] - step, myV, "myWall_Inner");
		verticalWall(comp, xSBmax[0] - step / 2,		lowY, yPad[0] - step, myV, "myWall_Inner");

		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - dXWall - step / 2,	topY, myV, "myWall_Inner");
		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - step,			lowY, myV, "myWall_Inner");

		verticalWall(comp, xMaxInnerGGWire - 0.5 - step, lowY, topY, myV, "myWall_Inner");

		comp->AddReadout("myWall_Inner");

//		Outer Wall
		verticalWall(comp, xSBmin[1] + dXWall - step / 2,	topY, yPad[1] - step / 2, myV, "myWall_Outer");
		verticalWall(comp, xSBmin[1] + step / 2,		lowY, yPad[1] - step / 2, myV, "myWall_Outer");

		horizontalWall(comp, xSBmin[1] + dXWall + step / 2,	xMinOuterGGWire + 0.5, topY, myV, "myWall_Outer");
		horizontalWall(comp, xSBmin[1] + 3 * step / 2,		xMinOuterGGWire + 0.5, lowY, myV, "myWall_Outer");

		verticalWall(comp, xMinOuterGGWire + 0.5 + step, lowY, topY, myV, "myWall_Outer");

		comp->AddReadout("myWall_Outer");
		}
	else if (Geometry.BeginsWith("myGeo_V02"))
		{
		double step = 0.011;

//		Inner Wall
		double myV = -190; 

//		position at myV
		double myY = (myV - vGG) / 133;

//		y positions of the walls
		double topY = myY + dXWall / 2 - step / 2;
		double lowY = myY - dXWall / 2 + step / 2;

		verticalWall(comp, xSBmax[0] - dXWall + step / 2,	topY, yPad[0] - step, myV, "myWall_Inner");
		verticalWall(comp, xSBmax[0] - step / 2,		lowY, yPad[0] - step, myV, "myWall_Inner");

		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - dXWall - step / 2,	topY, myV, "myWall_Inner");
		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - step,			lowY, myV, "myWall_Inner");

		verticalWall(comp, xMaxInnerGGWire - 0.5 - step, lowY, topY, myV, "myWall_Inner");

		comp->AddReadout("myWall_Inner");

//		Outer Wall
		verticalWall(comp, xSBmin[1] + dXWall - step / 2,	topY, yPad[1] - step / 2, myV, "myWall_Outer");
		verticalWall(comp, xSBmin[1] + step / 2,		lowY, yPad[1] - step / 2, myV, "myWall_Outer");

		horizontalWall(comp, xSBmin[1] + dXWall + step / 2,	xMinOuterGGWire + 0.5, topY, myV, "myWall_Outer");
		horizontalWall(comp, xSBmin[1] + 3 * step / 2,		xMinOuterGGWire + 0.5, lowY, myV, "myWall_Outer");

		verticalWall(comp, xMinOuterGGWire + 0.5 + step, lowY, topY, myV, "myWall_Outer");

		comp->AddReadout("myWall_Outer");
		}
	else if (Geometry.BeginsWith("myGeo_V03"))
		{
		dXWall = 0.08;
		double step = 0.011;

//		Inner Wall
		double myV = 0; 

//		position at myV
		double myY = (yFG + yGG) / 2;

//		y positions of the walls
		double topY = myY + dXWall / 2 - step / 2;
		double lowY = myY - dXWall / 2 + step / 2;

		verticalWall(comp, xSBmax[0] - dXWall + step / 2,	topY, yPad[0] - step, myV, "myWall_Inner");
		verticalWall(comp, xSBmax[0] - step / 2,		lowY, yPad[0] - step, myV, "myWall_Inner");

		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - dXWall - step / 2,	topY, myV, "myWall_Inner");
		horizontalWall(comp, xMaxInnerGGWire - 0.5, xSBmax[0] - step,			lowY, myV, "myWall_Inner");

		verticalWall(comp, xMaxInnerGGWire - 0.5 - step, lowY, topY, myV, "myWall_Inner");

		comp->AddReadout("myWall_Inner");

//		Outer Wall
		verticalWall(comp, xSBmin[1] + dXWall - step / 2,	topY, yPad[1] - step / 2, myV, "myWall_Outer");
		verticalWall(comp, xSBmin[1] + step / 2,		lowY, yPad[1] - step / 2, myV, "myWall_Outer");

		horizontalWall(comp, xSBmin[1] + dXWall + step / 2,	xMinOuterGGWire + 0.5, topY, myV, "myWall_Outer");
		horizontalWall(comp, xSBmin[1] + 3 * step / 2,		xMinOuterGGWire + 0.5, lowY, myV, "myWall_Outer");

		verticalWall(comp, xMinOuterGGWire + 0.5 + step, lowY, topY, myV, "myWall_Outer");

		comp->AddReadout("myWall_Outer");
		}
	else if (Geometry.BeginsWith("myGeo_V05"))
		{
		dXWall = 0.18;
		double step = 0.011;

//		position at myV
		double myY = (yFG + yGG) / 2;
		double myX = xGGmax[0] + 0.25; // 2.5 mm from the GG as requested by Jim;

//		wall outer side potential
		double outerV = -690;

		verticalWall(comp, myX + dXWall - step / 2,	myY, yPad[0] - step, outerV, "myWall");
		verticalWall(comp, myX + step / 2,		myY, yPad[0] - step, 0, "myWall");

		verticalWall(comp, myX + dXWall - step / 2,	yGG, myY - step, vGG, "myWall");
		verticalWall(comp, myX + step / 2,		yGG, myY - step, vGG, "myWall");

		horizontalWall(comp, myX + 1.5 * step, myX + dXWall - 1.5 * step, yGG, vGG, "myWall");

		comp->AddReadout("myWall");
		}
//	THIS IS A BRICK ON INNER INNER
	else if (Geometry.BeginsWith("myGeo_V06"))
		{
		dXWall = 0.18;
		double step = 0.011;

		double myX = xGGmin[0] + ggWpitch*2 - 0.3;
		horizontalWall(comp, myX - 0.6 + step, myX - step, yGG - step, 0, "myWall");
		verticalWall(comp, myX - 0.6, yGG, yPad[0] - step, 0, "myWall");
		verticalWall(comp, myX, yGG, yPad[0] - step, 0, "myWall");

		comp->AddReadout("myWall");
		}

//	======================= END OF MY GEOMETRY =========================================

//	Cathode plane
	comp->AddPlaneY(yC, vC, "q");
	comp->AddReadout("s");
	comp->AddReadout("p");
	comp->AddReadout("S");
	if (!  Geometry.BeginsWith("iTPC_3x125mkm_JT_092215_Proposal" ))
		{
		comp->AddReadout("Z");
		}
	comp->AddReadout("P");
//	comp->AddReadout("f");
//	comp->AddReadout("F");
	comp->SetGeometry(geo);
	TCanvas *c_e = 0, *c_i = 0;
	ViewDrift *v_e = 0, *v_i = 0;
	if (gROOT->IsBatch())
		{
		c_e = new TCanvas("Cell","Cell");
		TH1F *frame = c_e->DrawFrame(xmin-1.0,ymin-0.5,xmax+1.0,ymax+0.5);
		frame->SetXTitle("X (cm)");
		frame->SetYTitle("Y (cm)");
		frame->SetTitle("Cell structure");
		ViewCell *cellView = new ViewCell();
		cellView->SetComponent(comp);
		cellView->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
		cellView->SetCanvas(c_e);
		cellView->Plot2d();
		if (nEvents == 0) c_e->SaveAs("geometry.gif");
#if 0    
		c_i = new TCanvas("Cell_i","Cell for ions");
		ViewCell *cellView_i = new ViewCell();
		cellView_i->SetCanvas(c_i);
		cellView_i->SetComponent(comp);
		cellView_i->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
		cellView_i->Plot2d();
#endif
		}
//	Make a sensor
	Sensor* sensor = new Sensor();
//	sensor->EnableDebugging();
	sensor->AddComponent(comp);
	sensor->AddElectrode(comp, "p");
	sensor->SetTimeWindow(0., 50., 1000);
	sensor->ClearSignal();
	if (gROOT->IsBatch())
		{
		v_e = new ViewDrift();
		v_e->SetCanvas(c_e);
		v_e->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
		if (c_i)
			{
			v_i = new ViewDrift();
			v_i->SetCanvas(c_i);
			v_i->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
			}
//		Plot the potential
		ViewField* fView = new ViewField;
		fView->SetSensor(sensor);
		fView->SetArea(xmin, ymin, xmax, ymax);
		fView->SetVoltageRange(-200., 1400.);
//		fView->PlotSurface("e");
		TCanvas * fc = new TCanvas();
		fView->SetCanvas(fc);
		fView->PlotContour();
		if (nEvents == 0) fc->SaveAs("potential.gif");
		delete fc;
		}
	AvalancheMicroscopic* aval = new AvalancheMicroscopic(); 
//	aval->EnableDebugging();
	aval->EnableMagneticField();
	aval->SetSensor(sensor);
	if (! gROOT->IsBatch())
		{
		aval->EnablePlotting(v_e);
		}
//	aval->EnableNullCollisionSteps();
	DriftLineRKF* driftline_i = new DriftLineRKF();
	driftline_i->SetSensor(sensor);
	if (! gROOT->IsBatch())
		{
		driftline_i->EnablePlotting(v_i);
		}
	TH1::StatOverflows();
	TString outName;
	outName  = Region; outName += "_";
	outName += Geometry; outName += "_";
	outName += OutName;
	fOut = new TFile(outName,"recreate");
	if (! fOut) return;
	TH1F* hAngle = new TH1F("hAngle", "Angular distribution", nBinsAngular, -180, 180);
	TH2F* hRadial = new TH2F("hRadial", "Ion starting points", nBinsRadial, dSens / 2., 3. * dSens, nBinsAngular, -180, 180);
	TH1F* hElectrons = new TH1F("hElectrons", "Log_{10}(Number of electrons)", nBinsGain, 0, 5.);
	TH1F* hIons = new TH1F("hIons", "Log_{10}(Number of ions)", nBinsGain, 0, 5.);
	TH1F* tElectrons = new TH1F("tElectrons","time distribution of electors (ns)",100,0,1000);
	double step = 0.0100; 
	Int_t nx = (xmax - xmin)/step;
	Int_t ny = (ymax - ymin)/step;
	TH2F *xy_el   = new TH2F("el_xy"  ,"Electron at (x,y)", nx, xmin, xmax, ny, ymin, ymax);
	TH2F *xy_elL  = new TH2F("el_xyL" ,"Electron at (x,y) weighted with Length", nx, xmin, xmax, ny, ymin, ymax);
	TH2F *xy_elLT = new TH2F("el_xyLT","Electron at (x,y) weighted with Length and Time", nx, xmin, xmax, ny, ymin, ymax);
	TH2F *xy_ion   = new TH2F("ion_xy"  ,"Ion at (x,y)", nx, xmin, xmax, ny, ymin, ymax);
	TH2F *xy_ionL  = new TH2F("ion_xyL" ,"Ion at (x,y) weighted with Length", nx, xmin, xmax, ny, ymin, ymax);
	TH2F *xy_ionLT = new TH2F("ion_xyLT","Ion at (x,y) weighted with Length and Time", nx, xmin, xmax, ny, ymin, ymax);
	TH1F *el_status = new TH1F("el_status","electron status", 41, -20, 20);
	TH1F *ion_status = new TH1F("ion_status","ion status", 41, -20, 20);
	TH1F *eX = new TH1F("eX","Initial electrons",100,xmin,xmax);
	TProfile *eXI = new TProfile("eXI","Total no of ion versus x of primary electron",100,xmin,xmax);
	double x0, y0, z0, t0, e0;
	double x1, y1, z1, t1, e1;
	Int_t status;
	double r, phi;
	Int_t ne, ni;
	Int_t nEndpoints;
	for (Int_t i = nEvents; i--;)
		{
		if (i != nEvents) gBenchmark->Show("tpcGL");
		gBenchmark->Reset();
		gBenchmark->Start("tpcGL");
		gas->ResetCollisionCounters();
		sensor->NewSignal();

/*		THIS IS COMMENTED OUT AND BELOW BLOCK IS ENABLED TO REMOVE IONS ALONG WITH SEEDED ELECTRONS
		x0 = xmin + RndmUniform() * (xmax - xmin);
		eX->Fill(x0);
		y0 = yGG + 0.1;
		z0 = t0 = 0.;
		e0 = 0.1; // eV
		t0 = 0;
		aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);
*/
		double xStart = xmin + RndmUniform() * (xmax - xmin);
		eX->Fill(xStart);
		double yStart = yGG + 0.1;
		double zStart = 0;
		double tStart = 0.;
		double eStart = 0.1; // eV
		aval->AvalancheElectron(xStart, yStart, zStart, tStart, eStart, 0., 0., 0.);

//		v_e->Plot(true,false);
//		c_e->Update();
		aval->GetAvalancheSize(ne, ni);
		if (ne > 0) hElectrons->Fill(TMath::Log10(ne));    
		if (ni > 0) hIons->Fill(TMath::Log10(ni));
		eXI->Fill(x0,ni);
		nEndpoints = aval->GetNumberOfElectronEndpoints();
//		if (i % 10 == 0) 
		std::cout << i << "/" << nEvents << ": "
			<< ne << " electrons, "
			<< ni << " ions" << std::endl; 
		if (ne <= 0) continue;
		for (Int_t j = nEndpoints; j--;)
			{
			aval->GetElectronEndpoint(j, x0, y0, z0, t0, e0, x1, y1, z1, t1, e1, status);
//			BLOCK BELOW ADDED TO REMOVE IONS FROM SEEDED ELECTRONS -------
			double dist0 = TMath::Sqrt((x0 - xStart)*(x0 - xStart) + (y0 - yStart)*(y0 - yStart) + (z0 - zStart)*(z0 - zStart));
			double dist1 = TMath::Sqrt((x1 - xStart)*(x1 - xStart) + (y1 - yStart)*(y1 - yStart) + (z1 - zStart)*(z1 - zStart));
			if (dist0 < 1e-7 || dist1 < 1e-7) continue;
//			--------------------------------------------------------------

			tElectrons->Fill(t1);
			el_status->Fill(status);
			double X = (x1 + x0)/2;
			double Y = (y1 + y0)/2;
			double T = (t1 + t0)/2;
			double dX = x1 - x0;
			double dY = y1 - y0;
			double L = TMath::Sqrt(dX*dX + dY*dY);
//			if (j == 0) continue;
			xy_el->Fill(X,Y);
			xy_elL->Fill(X,Y,L);
			xy_elLT->Fill(X,Y,L*T);
			driftline_i->DriftIon(x0, y0, z0, t0);
			const std::vector<DriftLineRKF::step> &path_i = driftline_i->path();
			UInt_t nI = path_i.size();
			for (auto step : path_i)
				{
				ion_status->Fill(step.status);
				if (step.status) continue;
			        X = (step.xf + step.xi)/2;
			        Y = (step.yf + step.yi)/2;
			        T = (step.tf + step.ti)/2;
			        dX = step.xf - step.xi;
			        dY = step.yf - step.yi;
			        L = TMath::Sqrt(dX*dX + dY*dY);
				xy_ion->Fill(X,Y);
				xy_ionL->Fill(X,Y,L);
				xy_ionLT->Fill(X,Y,L*T);
				}
			}
//		v_i->Plot(true,false);
//		c_i->Update();
//		std::cout << "Next avalanche..." << std::endl;
		}
	gBenchmark->Show("tpcGL");
	fOut->Write();
	if (gROOT->IsBatch())  gSystem->Abort();
	}
//________________________________________________________________________________
void DrawFlux()
	{
	gStyle->SetOptStat(0);
	TH2F *ion_xyL = (TH2F *) gDirectory->Get("ion_xyL");
	if (! ion_xyL) return;
	double xmin = ion_xyL->GetXaxis()->GetXmin();
	double xmax = ion_xyL->GetXaxis()->GetXmax();
	TH1F *eX = (TH1F *) gDirectory->Get("eX");
	double N = 5000;
	if (eX)
		{
		N = eX->GetEntries();
		}
	double elFlux = N/(xmax - xmin);
	Int_t    nx   = ion_xyL->GetXaxis()->GetNbins();
	Int_t    ny   = ion_xyL->GetYaxis()->GetNbins();
	double ymin = ion_xyL->GetYaxis()->GetXmin();
	double ymax = ion_xyL->GetYaxis()->GetXmax();
	double scale = elFlux*(xmax - xmin)/nx*(ymax - ymin)/ny;
	cout << "scale = " << scale << endl;
	ion_xyL->Scale(1./scale);
	ion_xyL->SetTitle("Flux (ions/cm) normalized to electron flux 1 cm^{-1}");
	ion_xyL->SetXTitle("X (cm)");
	ion_xyL->SetYTitle("Z (cm)");
	ion_xyL->SetMinimum(1e-1);
	TCanvas *c1 = new TCanvas("fluxN","fluxN");
	c1->SetLogz(1);
	ion_xyL->Draw("colz");
	c1->SaveAs("fluxN.gif");
	TCanvas *c2 = new TCanvas("Projection","Projection");
	TH1D *proj = ion_xyL->ProjectionX("proj",11,40);
	proj->Scale(1./30);
	proj->SetTitle("Flux (ion/cm) at -0.25 cm from GG");
	proj->Draw();
	c2->SaveAs("Projection.gif");
	}
//________________________________________________________________________________
double stepF(double *x, double *p)
	{
	if (x[0] < p[1] - p[2]/2 || x[0] > p[1] + p[2]/2) return 0.;
	return p[0];
	}
//________________________________________________________________________________
TF1 *Step()
	{
	TF1 *f = new TF1("f",stepF,0,200,3);
	f->SetParameters(1.,100.,1.);
	f->SetParNames("IBF","X","FWHM");
	return f;
	}
/* 
      TPC:   IBF (ion/cm/e)      <X> (cm)  Width (cm)  
       OO:    8.100 +/- 0.300    192.00     1.0
1170 V II:    7.765 +/- 0.001     52.48     0.7
1100 V II:    5.320               52.49     0.7     
1170 V IO:    8.160              121.80     1.45
1100 V IO     7.223              121.80     1.45

OO/IO1170 = 8.100/8.160 = 0.99;   II1170/IO1170 = 7.765/8.160 = 0.95
OO/IO1110 = 8.100/7.223 = 1.12;   II1100/IO1100 = 5.320/7.223 = 0.74


 */
