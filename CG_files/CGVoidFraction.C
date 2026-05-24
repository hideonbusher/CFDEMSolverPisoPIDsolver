/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    This code is designed to realize coupled CFD-DEM simulations using LIGGGHTS
    and OpenFOAM(R). Note: this code is not part of OpenFOAM(R) (see DISCLAIMER).
\*---------------------------------------------------------------------------*/

#include "error.H"

#include "CGVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include <sstream>
#include <string>
#include "dictionary.H"
#include "IOobject.H"
#include "List.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CGVoidFraction, 0);

addToRunTimeSelectionTable
(
    voidFractionModel,
    CGVoidFraction,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
CGVoidFraction::CGVoidFraction
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    voidFractionModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
	mesh_(particleCloud_.mesh()),
    verbose_(false),
    procBoundaryCorrection_(propsDict_.lookupOrDefault<Switch>("procBoundaryCorrection", false)),
    alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
	meshx_(readScalar(propsDict_.lookup("meshx"))),
	meshy_(readScalar(propsDict_.lookup("meshy"))),
	meshz_(readScalar(propsDict_.lookup("meshz"))),
	decx_(readScalar(propsDict_.lookup("decx"))),
	decy_(readScalar(propsDict_.lookup("decy"))),
	decz_(readScalar(propsDict_.lookup("decz"))),
	periodicx_(propsDict_.lookupOrDefault<Switch>("periodicx", false)),
    periodicy_(propsDict_.lookupOrDefault<Switch>("periodicy", false)),
    periodicz_(propsDict_.lookupOrDefault<Switch>("periodicz", false)),
	TopologyAccel_(propsDict_.lookupOrDefault<Switch>("TopologyAccel", false)),
    alphaLimited_(0),
    tooMuch_(0.0),
    interpolation_(false),
    cfdemUseOnly_(false)
{
    maxCellsPerParticle_ = 313;
    //particleCloud_.setMaxCellsPerParticle(29);

    if(alphaMin_ > 1 || alphaMin_ < 0.01){ FatalError<< "alphaMin should be < 1 and > 0.01 !!!" << abort(FatalError); }
    if (propsDict_.found("interpolation")){
        interpolation_=true;
        Warning << "interpolation for CGVoidFraction does not yet work correctly!" << endl;
        //Info << "Using interpolated voidfraction field - do not use this in combination with interpolation in drag model!"<< endl;
    }

    checkWeightNporosity(propsDict_);

    if (propsDict_.found("verbose")) verbose_=true;

    if (propsDict_.found("cfdemUseOnly"))
    {
        cfdemUseOnly_ = readBool(propsDict_.lookup("cfdemUseOnly"));
    }

    // check if settings are consistent with locate model selected
    if (procBoundaryCorrection_)
    {
        if(!(particleCloud_.locateM().type()=="engineIB"))
        {
            FatalError << typeName << ": You are requesting procBoundaryCorrection, this requires the use of engineIB!\n"
                       << abort(FatalError);
        }
    } else {
        if(particleCloud_.locateM().type()=="engineIB")
        {
            FatalError << typeName << ": You are using engineIB, this requires using procBoundaryCorrection=true!\n"
                       << abort(FatalError);
            //Warning << "You are trying to use engineIB, this requires using procBoundaryCorrection=true\n"
            //        << "  procBoundaryCorrection will be used!\n" << endl;
            //procBoundaryCorrection_ = true;
        }//
    }
	
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CGVoidFraction::~CGVoidFraction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label CGVoidFraction::calculateLocalGridID(
    const vector& subPosition,
    const boundBox& globalBb,
    int meshx, int meshy, int meshz,
    scalar decx_, scalar decy_, scalar decz_) const
{
    // Compute grid spacing
    scalar dx = (globalBb.max()[0] - globalBb.min()[0]) / meshx; // grid spacing in the x direction
    scalar dy = (globalBb.max()[1] - globalBb.min()[1]) / meshy; // grid spacing in the y direction
    scalar dz = (globalBb.max()[2] - globalBb.min()[2]) / meshz; // grid spacing in the z direction

    // Convert decx_, decy_, and decz_ to integers
    int decx = static_cast<int>(decx_);
    int decy = static_cast<int>(decy_);
    int decz = static_cast<int>(decz_);

    // Current processor ID and its indices in the x, y, and z directions
    int procID = Pstream::myProcNo();
    int procX = procID % decx;
    int procY = (procID / decx) % decy;
    int procZ = procID / (decx * decy);

    // Number of grid cells handled by the local processor
    int ProMeshNumX = meshx / decx;
    int ProMeshNumY = meshy / decy;
    int ProMeshNumZ = meshz / decz;

    // Local processor range
    scalar localProcessorXMin = globalBb.min()[0] + procX * ProMeshNumX * dx;
    scalar localProcessorXMax = localProcessorXMin + ProMeshNumX * dx;

    scalar localProcessorYMin = globalBb.min()[1] + procY * ProMeshNumY * dy;
    scalar localProcessorYMax = localProcessorYMin + ProMeshNumY * dy;

    scalar localProcessorZMin = globalBb.min()[2] + procZ * ProMeshNumZ * dz;
    scalar localProcessorZMax = localProcessorZMin + ProMeshNumZ * dz;

    // Check whether the point is within the current processor range
    if (subPosition[0] < localProcessorXMin || subPosition[0] >= localProcessorXMax ||
        subPosition[1] < localProcessorYMin || subPosition[1] >= localProcessorYMax ||
        subPosition[2] < localProcessorZMin || subPosition[2] >= localProcessorZMax)
    {
        return -1; // The point is not within the current processor range
    }

    // Compute local grid indices
    int I = static_cast<int>((subPosition[0] - localProcessorXMin) / dx);
    int J = static_cast<int>((subPosition[1] - localProcessorYMin) / dy);
    int K = static_cast<int>((subPosition[2] - localProcessorZMin) / dz);

    // Compute the local grid ID according to the index formula
    //return K * (ProMeshNumX * ProMeshNumY) + (ProMeshNumX - 1 - I) * ProMeshNumY + J;
	return K * (ProMeshNumX * ProMeshNumY) + J * ProMeshNumX + I;
}



void CGVoidFraction::setvoidFraction(double** const& mask,double**& voidfractions,double**& particleWeights,double**& particleVolumes, double**& particleV) const
{
    if(cfdemUseOnly_)
        reAllocArrays(particleCloud_.numberOfParticles());
    else
        reAllocArrays();
    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);//0826
	//for(int aa=0; aa< 1000; aa++)
	//{Info << " aa=" << aa << endl; Info << " voidfractionNext_=" << voidfractionNext_[aa] << endl; }
	
    vector position(0.,0.,0.);
    label cellID = -1;
    scalar radius(-1.);
    scalar volume(0.);
    scalar cellVol(0.);
    scalar scaleVol= weight();
    scalar scaleRadius = cbrt(porosity());
	scalar ds(0);
    const boundBox& globalBb = particleCloud_.mesh().bounds();
    label partCellId = -1 ;	
////////////////////////////////////Set satellite-point coefficients//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	scalar r[] = { 0.5,1,1.5,2,2.5,3,3.5,4 };//Define satellite-point layers 
	vector offsets[numberOfMarkerPoints];// add
    int m = 0;
    offsets[m][0] = offsets[m][1] = offsets[m][2] = 0.0;//Initialize and assign the center point
    m ++;//m = satellite point	
	double phipi=(sqrt(5.)-1)/2;
	    //0.5 layer: 6 points
        for (int j = 1; j <= 6; j++)
	    {
			offsets[m][2] = r[0]*((2*static_cast<double>(j)-1)/6-1);
			offsets[m][0] = r[0]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[0]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
        //First layer: 6 points
        for (int j = 1; j <= 6; j++)
	    {
			offsets[m][2] = r[1]*((2*static_cast<double>(j)-1)/6-1);
			offsets[m][0] = r[1]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[1]*sqrt(1-((2*static_cast<double>(j)-1)/6-1)*((2*static_cast<double>(j)-1)/6-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
		//1.5 layer: 14 points
        for (int j = 1; j <= 14; j++)
	    {
			offsets[m][2] = r[2]*((2*static_cast<double>(j)-1)/14-1);
			offsets[m][0] = r[2]*sqrt(1-((2*static_cast<double>(j)-1)/14-1)*((2*static_cast<double>(j)-1)/14-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[2]*sqrt(1-((2*static_cast<double>(j)-1)/14-1)*((2*static_cast<double>(j)-1)/14-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
			m++;
		}
		//Second layer: 24 points
        for (int j = 1; j <= 24; j++) 
        {
			offsets[m][2] = r[3]*((2*static_cast<double>(j)-1)/24-1);
			offsets[m][0] = r[3]*sqrt(1-((2*static_cast<double>(j)-1)/24-1)*((2*static_cast<double>(j)-1)/24-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[3]*sqrt(1-((2*static_cast<double>(j)-1)/24-1)*((2*static_cast<double>(j)-1)/24-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m2=" << m << endl;
        }
		//2.5 layer: 38 points
        for (int j = 1; j <= 38; j++) 
        {
			offsets[m][2] = r[4]*((2*static_cast<double>(j)-1)/38-1);
			offsets[m][0] = r[4]*sqrt(1-((2*static_cast<double>(j)-1)/38-1)*((2*static_cast<double>(j)-1)/38-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[4]*sqrt(1-((2*static_cast<double>(j)-1)/38-1)*((2*static_cast<double>(j)-1)/38-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m2=" << m << endl;
        }
		//Third layer: 54 points
		for (int j = 1; j <= 54; j++) 
        {
			offsets[m][2] = r[5]*((2*static_cast<double>(j)-1)/54-1);
			offsets[m][0] = r[5]*sqrt(1-((2*static_cast<double>(j)-1)/54-1)*((2*static_cast<double>(j)-1)/54-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[5]*sqrt(1-((2*static_cast<double>(j)-1)/54-1)*((2*static_cast<double>(j)-1)/54-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m3=" << m << endl;
        }
		//3.5 layer: 74 points
		for (int j = 1; j <= 74; j++) 
        {
			offsets[m][2] = r[6]*((2*static_cast<double>(j)-1)/74-1);
			offsets[m][0] = r[6]*sqrt(1-((2*static_cast<double>(j)-1)/74-1)*((2*static_cast<double>(j)-1)/74-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[6]*sqrt(1-((2*static_cast<double>(j)-1)/74-1)*((2*static_cast<double>(j)-1)/74-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m3=" << m << endl;
        }
		//Fourth layer: 96 points
		for (int j = 1; j <= 96; j++) 
        {
			offsets[m][2] = r[7]*((2*static_cast<double>(j)-1)/96-1);
			offsets[m][0] = r[7]*sqrt(1-((2*static_cast<double>(j)-1)/96-1)*((2*static_cast<double>(j)-1)/96-1))*Foam::cos(2*M_PI*phipi*static_cast<double>(j));
            offsets[m][1] = r[7]*sqrt(1-((2*static_cast<double>(j)-1)/96-1)*((2*static_cast<double>(j)-1)/96-1))*Foam::sin(2*M_PI*phipi*static_cast<double>(j));
            m ++;//Pout << " m4=" << m << endl;
        }
//#pragma omp parallel for
if (TopologyAccel_)
{
	for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {

        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        //if(mask[index][0])
        //{
            // reset

            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0.;
                particleVolumes[index][subcell] = 0.;
            }
            particleV[index][0] = 0.;

            cellsPerParticle_[index][0] = 1.;//cellsPerParticle is the number of cells occupied (overlapped) by the particle
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
			radius = particleRadius(index);//particleCloud_.radius(index);
            volume = Vp(index,radius,scaleVol);
            radius *= scaleRadius;
            cellVol = 0.;
            ds = 2*particleCloud_.radius(index);
            //--variables for sub-search
            int nPoints = numberOfMarkerPoints; //numberOfMarkerPoints is defined in the header file
            int nNotFound=0,nUnEqual=0,nTotal=0;
            vector offset(0.,0.,0.);
            int cellsSet = 0;
            label cellWithCenter(-1);
			//In OpenFOAM, label is a specific integer type used to identify mesh elements such as cells, faces, edges, and points
            //The label type in OpenFOAM is usually defined as a 32-bit integer (int), although this may depend on the version or compilation options
			//Variables of type label are special integer identifiers for mesh elements in OpenFOAM, with the above characteristics and uses
			if(procBoundaryCorrection_)
            {
                // switch off cellIDs for force calc if steming from parallel search success
                cellWithCenter = particleCloud_.locateM().findSingleCell(position,cellID);
                particleCloud_.cellIDs()[index][0] = cellWithCenter;
            }
	labelList transffi(313,0);
	int count=0;
	int cc=0;
/*///////////////Local-global cell ID/////////////////////////////////////////////////////
#include "globalMeshData.H"
#include "globalIndex.H"
const polyMesh& mesh = particleCloud_.mesh();
Foam::globalMeshData globalData(mesh); // Create a globalMeshData object
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering(); // Get the global cell-numbering utility
label globalID = globalCellIndex.toGlobal(cellID);
// Output the result
Pout << "Processor " << Foam::Pstream::myProcNo()
     << ", Local Cell ID: " << cellID
     << ", Global Cell ID: " << globalID << Foam::endl;
	 	
#include "Pstream.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "HashSet.H"

const polyMesh& mesh = particleCloud_.mesh();
const Foam::globalMeshData globalData(mesh);

// Get the global-numbering utility
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering();

// Get offsets
const Foam::labelList& ooffsets = const_cast<Foam::globalIndex&>(globalCellIndex).offsets();

Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
           << ", Offsets: " << ooffsets << Foam::endl;


/*for (Foam::label localID = 0; localID < mesh.nCells(); ++localID)
{
    Foam::label globalID = globalCellIndex.toGlobal(localID);

    // Determine whether the current processor is the master processor
    if (globalCellIndex.whichProcID(globalID) == Foam::Pstream::myProcNo())
    {
        Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
                   << ", Local Cell ID: " << localID
                   << ", Global Cell ID: " << globalID << " (Master)" << Foam::endl;
    }
}*/


/*const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
forAll(boundaryMesh, patchi)
{
    const polyPatch& patch = boundaryMesh[patchi];

    if (isA<processorPolyPatch>(patch))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
        const labelList& faceCells = procPatch.faceCells();

        // Loop over shared cells
        forAll(faceCells, i)
        {
            Foam::Pout << "Shared Cell ID: " << faceCells[i] << Foam::endl;
        }
    }
}*/

//////////////////////////////////////////////////////////////////
 if (cellID >= 0)  
 {
    cellVol = particleCloud_.mesh().V()[cellID];
	scalar lostweight=0;
  for(int i = 0; i < numberOfMarkerPoints; i++) //Loop over all satellite points and compute their offsets
  {

//////////////////////////////////////////CGprocess//////////////////////////////////////////////////////////////////						                 
    // locate subPoint
	                vector offset(0.,0.,0.);
					offset = particleRadius(index)*offsets[i];
					scalar distance = mag (offset);
					vector subPosition = particleCloud_.position(index) + offset;
					cellVol = particleCloud_.mesh().V()[cellID];
					scalar wl = 0;
					scalar layerR=0;

		// Apply periodic correction to each direction of subPosition, iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // Check whether periodic boundary conditions are enabled in the current direction
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}

	  
					
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
	//if (TopologyAccel_)			
    label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
    //else
	//{label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
	//label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
	//label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
    //Pout << "partCellId true =" << partCellIdtrue << " ----partCellId TST = " << partCellId <<  endl;
     if (verbose_ && index==0)  	 
		 {
				label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                Pout << "partCellId true =" << partCellIdtrue << "--V-S--partCellId TST = " << partCellId <<  endl;
		 }	
////////////////////////////virtual subparticle method liuyuxiang2024.9.16	
    if (partCellId >= 0)  // subPoint is in domain check whether the satellite point is inside the computational domain
    {
        // update voidfraction for each particle read
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];//extract the volume of the cell containing this satellite point
        scalar particleVolume = volume*wl;//weighted volume of the satellite point
      
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++;   
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
        bool createNew = true;
        label storeInIndex=0;
        for(int i=0; i < cellsPerParticle_[index][0] ; i++)
        {
            if(partCellId == particleCloud_.cellIDs()[index][i]) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }
        }
        if(createNew)
        {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0]-1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
        }
        particleWeights[index][storeInIndex] += wl; 
        particleVolumes[index][storeInIndex] += particleVolume;  
        particleV[index][0] += particleVolume;	
    }
	if (partCellId < 0) 
	{    
        ////////////////////////////////////////first physical boundary
		if (subPosition[0] > globalBb.max()[0] || subPosition[1] > globalBb.max()[1] || subPosition[2] > globalBb.max()[2] || subPosition[0] < globalBb.min()[0] || subPosition[1] < globalBb.min()[1] || subPosition[2] < globalBb.min()[2])
	    {    
            while (partCellId < 0)//&& iterations < iterationslimit)	//virtual sub-point
            { 
				if (layerR <= 0) //treatment for point inside processor boundary!	(very little)
				{
					break;
				}		
	            for (int ivsp = 0 ; ivsp < 3 ; ivsp++)
	            {
                    offset[ivsp] = offset[ivsp]*(layerR-0.5)/(layerR);// layerR must not be 0
	            }	
	            layerR=layerR-0.5;
	            subPosition = particleCloud_.position(index) + offset;	
	            //partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
				//partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
    			//if (TopologyAccel_)			
   			    partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
   			    //else
				//{partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
            }
			
			if (partCellId >= 0) //////// move the virtual point into computational domain
			{   
			    scalar partCellVol = particleCloud_.mesh().V()[partCellId];//extract the volume of the cell containing this satellite point
                scalar particleVolume = volume*wl;//weighted volume of the satellite point
                scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		        
        	    if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        	    else
        	    {
        	        voidfractionNext_[partCellId] = alphaMin_;
        	        tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        	    }
        	    cellsSet++; 
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
				bool createNew = true;
 				label storeInIndex=0;
  				for(int i=0; i < cellsPerParticle_[index][0] ; i++)
  				{
    				if(partCellId == particleCloud_.cellIDs()[index][i]) 
    				{
      				    storeInIndex = i;
     				    createNew = false;
      				    break;
     				}
     			}
     			if(createNew)
     			{
      				cellsPerParticle_[index][0] ++;
      				storeInIndex = cellsPerParticle_[index][0]-1;
      				particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
      			}
      			particleWeights[index][storeInIndex] += wl; 
      			particleVolumes[index][storeInIndex] += particleVolume;  
       			particleV[index][0] += particleVolume;
			}
	    }//end physical boundary revision (out of physcial boundary!)
		else // inside the physcial boundary but not on this processor!
		{ 
		    //Weight correction: the void fraction is accurately transferred by reduce, while the weights are corrected only on this processor
			lostweight= lostweight+wl; //Correct only the weights, not the void fraction, because the void fraction is handled by reduce later
		    transffi[cc]=i;
		    cc=cc+1;
		    count=count+1;
		}
    }//end for not main processor (partcellid not found!)
 }// end subpoint loop
    for(int subcellss=0;subcellss<cellsPerParticle_[index][0];subcellss++)
    {
	  particleWeights[index][subcellss]=particleWeights[index][subcellss]/(1-lostweight);
	  particleVolumes[index][subcellss]=particleVolumes[index][subcellss]/(1-lostweight);
    }
	particleV[index][0]=particleV[index][0]/(1-lostweight);
}//end if in cell main processor

	 reduce(transffi, maxOp<List<label>>());// test OK transffi=i corresponding to the unfound partID (point-cloud order)
	 reduce(count, maxOp<int>());// test OK transffi=i corresponding to the unfound partID (point-cloud order)
if (cellID < 0)
{
	for (cc=0;cc<count;cc++) //Loop over satellite points not on the main processor on non-main processors
	{
		int i=transffi[cc]; // transffi stores i in ascending order; the point-cloud order is equivalent to i on the main processor, but here i corresponds to points that were not found
		vector offset(0.,0.,0.);
		offset = particleRadius(index)*offsets[i];
		scalar distance = mag (offset);
		vector subPosition = particleCloud_.position(index) + offset;
		cellVol = particleCloud_.mesh().V()[cellID];
		scalar wl = 0;
		scalar layerR=0;
// Apply periodic correction to each direction of subPosition, iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // Check whether periodic boundary conditions are enabled in the current direction
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
 		   // If this direction is periodic, wrap the coordinate
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
                    //label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                    //label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
				//if (TopologyAccel_)			
   			    label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); 
   			    //else
				//{label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }
	  if (partCellId >= 0)///////////////satellite points on non-main processors
	  {
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];
        scalar particleVolume = volume*wl;
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++; //cancel weighting! because the error of array only revise the voidfraction!
	  } 
	}// end subpoint loop for satellite points on non-main processors		 
}// end non-main processor	

	            if (verbose_ && index == 0)
				{
				   scalar sum=0; 
				   for (int lyx = 0; lyx < particleCloud_.cellsPerParticle()[index][0]; lyx++) //kernelcell is the cell sequence occupied by the particle CG domain, while storesubcellkernal is the actual cell ID
				   {
				       sum=sum+particleCloud_.particleWeights()[index][lyx]; //debug for mass conservation
                   }
				   Pout << "volume=" << volume << endl; 
				   Pout << "particleWeights[index][0]=" << particleWeights[index][0] << endl;
				   Pout << "particleVolumes[index][0]=" << particleVolumes[index][0] << endl;
				   Pout << "particleV[index][0]=" << particleV[index][0] << endl;
				   Pout << "sum=" << sum << endl;
				}

    }//// end loop all particles
}//// end TopologyAccel_
else
{		
    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {

        if(!checkParticleType(index)) continue; //skip this particle if not correct type

        //if(mask[index][0])
        //{
            // reset

            for(int subcell=0;subcell<cellsPerParticle_[index][0];subcell++)
            {
                particleWeights[index][subcell] = 0.;
                particleVolumes[index][subcell] = 0.;
            }
            particleV[index][0] = 0.;

            cellsPerParticle_[index][0] = 1.;//cellsPerParticle is the number of cells occupied (overlapped) by the particle
            position = particleCloud_.position(index);
            cellID = particleCloud_.cellIDs()[index][0];
			radius = particleRadius(index);//particleCloud_.radius(index);
            volume = Vp(index,radius,scaleVol);
            radius *= scaleRadius;
            cellVol = 0.;
            ds = 2*particleCloud_.radius(index);
            //--variables for sub-search
            int nPoints = numberOfMarkerPoints; //numberOfMarkerPoints is defined in the header file
            int nNotFound=0,nUnEqual=0,nTotal=0;
            vector offset(0.,0.,0.);
            int cellsSet = 0;
            label cellWithCenter(-1);
			//In OpenFOAM, label is a specific integer type used to identify mesh elements such as cells, faces, edges, and points
            //The label type in OpenFOAM is usually defined as a 32-bit integer (int), although this may depend on the version or compilation options
			//Variables of type label are special integer identifiers for mesh elements in OpenFOAM, with the above characteristics and uses
			if(procBoundaryCorrection_)
            {
                // switch off cellIDs for force calc if steming from parallel search success
                cellWithCenter = particleCloud_.locateM().findSingleCell(position,cellID);
                particleCloud_.cellIDs()[index][0] = cellWithCenter;
            }
	labelList transffi(313,0);
	int count=0;
	int cc=0;
/*///////////////Local-global cell ID/////////////////////////////////////////////////////
#include "globalMeshData.H"
#include "globalIndex.H"
const polyMesh& mesh = particleCloud_.mesh();
Foam::globalMeshData globalData(mesh); // Create a globalMeshData object
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering(); // Get the global cell-numbering utility
label globalID = globalCellIndex.toGlobal(cellID);
// Output the result
Pout << "Processor " << Foam::Pstream::myProcNo()
     << ", Local Cell ID: " << cellID
     << ", Global Cell ID: " << globalID << Foam::endl;
	 	
#include "Pstream.H"
#include "globalMeshData.H"
#include "globalIndex.H"
#include "HashSet.H"

const polyMesh& mesh = particleCloud_.mesh();
const Foam::globalMeshData globalData(mesh);

// Get the global-numbering utility
const Foam::globalIndex& globalCellIndex = globalData.globalBoundaryCellNumbering();

// Get offsets
const Foam::labelList& ooffsets = const_cast<Foam::globalIndex&>(globalCellIndex).offsets();

Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
           << ", Offsets: " << ooffsets << Foam::endl;


/*for (Foam::label localID = 0; localID < mesh.nCells(); ++localID)
{
    Foam::label globalID = globalCellIndex.toGlobal(localID);

    // Determine whether the current processor is the master processor
    if (globalCellIndex.whichProcID(globalID) == Foam::Pstream::myProcNo())
    {
        Foam::Pout << "Processor " << Foam::Pstream::myProcNo()
                   << ", Local Cell ID: " << localID
                   << ", Global Cell ID: " << globalID << " (Master)" << Foam::endl;
    }
}*/


/*const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
forAll(boundaryMesh, patchi)
{
    const polyPatch& patch = boundaryMesh[patchi];

    if (isA<processorPolyPatch>(patch))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
        const labelList& faceCells = procPatch.faceCells();

        // Loop over shared cells
        forAll(faceCells, i)
        {
            Foam::Pout << "Shared Cell ID: " << faceCells[i] << Foam::endl;
        }
    }
}*/

//////////////////////////////////////////////////////////////////
 if (cellID >= 0)  
 {
    cellVol = particleCloud_.mesh().V()[cellID];
	scalar lostweight=0;
  for(int i = 0; i < numberOfMarkerPoints; i++) //Loop over all satellite points and compute their offsets
  {

//////////////////////////////////////////CGprocess//////////////////////////////////////////////////////////////////						                 
    // locate subPoint
	                vector offset(0.,0.,0.);
					offset = particleRadius(index)*offsets[i];
					scalar distance = mag (offset);
					vector subPosition = particleCloud_.position(index) + offset;
					cellVol = particleCloud_.mesh().V()[cellID];
					scalar wl = 0;
					scalar layerR=0;

		// Apply periodic correction to each direction of subPosition, iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // Check whether periodic boundary conditions are enabled in the current direction
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}

	  
					
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
	//if (TopologyAccel_)			
    //{label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
    //else
	label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); 	
	//label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
	//label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
    //Pout << "partCellId true =" << partCellIdtrue << " ----partCellId TST = " << partCellId <<  endl;
     if (verbose_ && index==0)  	 
		 {
				label partCellIdtrue = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                Pout << "partCellId true =" << partCellIdtrue << "--V-S--partCellId TST = " << partCellId <<  endl;
		 }	
////////////////////////////virtual subparticle method liuyuxiang2024.9.16	
    if (partCellId >= 0)  // subPoint is in domain check whether the satellite point is inside the computational domain
    {
        // update voidfraction for each particle read
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];//extract the volume of the cell containing this satellite point
        scalar particleVolume = volume*wl;//weighted volume of the satellite point
      
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++;   
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
        bool createNew = true;
        label storeInIndex=0;
        for(int i=0; i < cellsPerParticle_[index][0] ; i++)
        {
            if(partCellId == particleCloud_.cellIDs()[index][i]) 
            {
                storeInIndex = i;
                createNew = false;
                break;
            }
        }
        if(createNew)
        {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0]-1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
        }
        particleWeights[index][storeInIndex] += wl; 
        particleVolumes[index][storeInIndex] += particleVolume;  
        particleV[index][0] += particleVolume;	
    }
	if (partCellId < 0) 
	{    
        ////////////////////////////////////////first physical boundary
		if (subPosition[0] > globalBb.max()[0] || subPosition[1] > globalBb.max()[1] || subPosition[2] > globalBb.max()[2] || subPosition[0] < globalBb.min()[0] || subPosition[1] < globalBb.min()[1] || subPosition[2] < globalBb.min()[2])
	    {    
            while (partCellId < 0)//&& iterations < iterationslimit)	//virtual sub-point
            { 
				if (layerR <= 0) //treatment for point inside processor boundary!	(very little)
				{
					break;
				}		
	            for (int ivsp = 0 ; ivsp < 3 ; ivsp++)
	            {
                    offset[ivsp] = offset[ivsp]*(layerR-0.5)/(layerR);// layerR must not be 0
	            }	
	            layerR=layerR-0.5;
	            subPosition = particleCloud_.position(index) + offset;	
	            //partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
				//partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
    			//if (TopologyAccel_)			
   			    //{partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
   			    //else
				{partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); }	
            }
			
			if (partCellId >= 0) //////// move the virtual point into computational domain
			{   
			    scalar partCellVol = particleCloud_.mesh().V()[partCellId];//extract the volume of the cell containing this satellite point
                scalar particleVolume = volume*wl;//weighted volume of the satellite point
                scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
		        
        	    if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        	    else
        	    {
        	        voidfractionNext_[partCellId] = alphaMin_;
        	        tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        	    }
        	    cellsSet++; 
//--------set sub weight--------////--------set sub weight--------////--------set sub weight--------////--------set sub weight--------//
				bool createNew = true;
 				label storeInIndex=0;
  				for(int i=0; i < cellsPerParticle_[index][0] ; i++)
  				{
    				if(partCellId == particleCloud_.cellIDs()[index][i]) 
    				{
      				    storeInIndex = i;
     				    createNew = false;
      				    break;
     				}
     			}
     			if(createNew)
     			{
      				cellsPerParticle_[index][0] ++;
      				storeInIndex = cellsPerParticle_[index][0]-1;
      				particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
      			}
      			particleWeights[index][storeInIndex] += wl; 
      			particleVolumes[index][storeInIndex] += particleVolume;  
       			particleV[index][0] += particleVolume;
			}
	    }//end physical boundary revision (out of physcial boundary!)
		else // inside the physcial boundary but not on this processor!
		{ 
		    //Weight correction: the void fraction is accurately transferred by reduce, while the weights are corrected only on this processor
			lostweight= lostweight+wl; //Correct only the weights, not the void fraction, because the void fraction is handled by reduce later
		    transffi[cc]=i;
		    cc=cc+1;
		    count=count+1;
		}
    }//end for not main processor (partcellid not found!)
 }// end subpoint loop
    for(int subcellss=0;subcellss<cellsPerParticle_[index][0];subcellss++)
    {
	  particleWeights[index][subcellss]=particleWeights[index][subcellss]/(1-lostweight);
	  particleVolumes[index][subcellss]=particleVolumes[index][subcellss]/(1-lostweight);
    }
	particleV[index][0]=particleV[index][0]/(1-lostweight);
}//end if in cell main processor

	 reduce(transffi, maxOp<List<label>>());// test OK transffi=i corresponding to the unfound partID (point-cloud order)
	 reduce(count, maxOp<int>());// test OK transffi=i corresponding to the unfound partID (point-cloud order)
if (cellID < 0)
{
	for (cc=0;cc<count;cc++) //Loop over satellite points not on the main processor on non-main processors
	{
		int i=transffi[cc]; // transffi stores i in ascending order; the point-cloud order is equivalent to i on the main processor, but here i corresponds to points that were not found
		vector offset(0.,0.,0.);
		offset = particleRadius(index)*offsets[i];
		scalar distance = mag (offset);
		vector subPosition = particleCloud_.position(index) + offset;
		cellVol = particleCloud_.mesh().V()[cellID];
		scalar wl = 0;
		scalar layerR=0;
// Apply periodic correction to each direction of subPosition, iDir: 0->x, 1->y, 2->z
		for (int iDir = 0; iDir < 3; iDir++)
		{
 		   // Check whether periodic boundary conditions are enabled in the current direction
 		   bool isPeriodic = false;
 		   if (iDir == 0 && periodicx_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 1 && periodicy_)
 		   {
 		       isPeriodic = true;
 		   }
 		   else if (iDir == 2 && periodicz_)
 		   {
  		      isPeriodic = true;
 		   }
		
 		   // If this direction is periodic, wrap the coordinate
  		    if (isPeriodic)
  		   {
  		      if (subPosition[iDir] > globalBb.max()[iDir])
   		     {
           		 subPosition[iDir] -= (globalBb.max()[iDir] - globalBb.min()[iDir]);
        		}
 		       else if (subPosition[iDir] < globalBb.min()[iDir])
   		     {
   		         subPosition[iDir] += (globalBb.max()[iDir] - globalBb.min()[iDir]);
  		      }
   		    }
		}
                        if (i == 0)
						{ wl=0.130856838989584*exp(-distance*distance/(8*ds*ds)); }
					
					else if (i >= 1 && i <= 6)
					{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 0.5 ;}
				
					    else if (i >= 7 && i <= 12)
						{ wl=0.130856838989584/6*exp(-distance*distance/(8*ds*ds)); layerR = 1 ;}
					
					else if (i >= 13 && i <= 26)
					{ wl=0.130856838989584/14*exp(-distance*distance/(8*ds*ds)); layerR = 1.5 ;}
				
					    else if (i >= 27 && i <= 50)
						{ wl=0.130856838989584/24*exp(-distance*distance/(8*ds*ds)); layerR = 2 ;}
					
					else if (i >= 51 && i <= 88)
					{ wl=0.130856838989584/38*exp(-distance*distance/(8*ds*ds)); layerR = 2.5 ;}
				
					    else if (i >= 89 && i <= 142)
						{ wl=0.130856838989584/54*exp(-distance*distance/(8*ds*ds)); layerR = 3 ;}
					
					else if (i >= 143 && i <= 216)
					{ wl=0.130856838989584/74*exp(-distance*distance/(8*ds*ds)); layerR = 3.5 ;}
				
					    else if (i >= 217 && i <= 312)
						{ wl=0.130856838989584/96*exp(-distance*distance/(8*ds*ds)); layerR = 4 ;}
                    //label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID);
                    //label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_);
				//if (TopologyAccel_)			
   			    //{label partCellId = calculateLocalGridID(subPosition, globalBb, meshx_, meshy_, meshz_, decx_, decy_, decz_); }
   			    //else
				label partCellId = particleCloud_.locateM().findSingleCell(subPosition,cellID); 
	  if (partCellId >= 0)///////////////satellite points on non-main processors
	  {
        scalar partCellVol = particleCloud_.mesh().V()[partCellId];
        scalar particleVolume = volume*wl;
        scalar newAlpha = voidfractionNext_[partCellId]- particleVolume / partCellVol;
        if(newAlpha > alphaMin_) voidfractionNext_[partCellId] = newAlpha;
        else
        {
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_-newAlpha) * partCellVol;
        }
        cellsSet++; //cancel weighting! because the error of array only revise the voidfraction!
	  } 
	}// end subpoint loop for satellite points on non-main processors		 
}// end non-main processor	

	            if (verbose_ && index == 0)
				{
				   scalar sum=0; 
				   for (int lyx = 0; lyx < particleCloud_.cellsPerParticle()[index][0]; lyx++) //kernelcell is the cell sequence occupied by the particle CG domain, while storesubcellkernal is the actual cell ID
				   {
				       sum=sum+particleCloud_.particleWeights()[index][lyx]; //debug for mass conservation
                   }
				   Pout << "volume=" << volume << endl; 
				   Pout << "particleWeights[index][0]=" << particleWeights[index][0] << endl;
				   Pout << "particleVolumes[index][0]=" << particleVolumes[index][0] << endl;
				   Pout << "particleV[index][0]=" << particleV[index][0] << endl;
				   Pout << "sum=" << sum << endl;
				}

    }//// end loop all particles
}





///////////////new loop test/////////////////////
    voidfractionNext_.correctBoundaryConditions();
	//for(int aa=0; aa< 1000; aa++)
	//{Info << voidfractionNext_[aa] << " voidfraction=" << endl; }	
    //Used to correct or ensure that the values in voidfractionNext_ satisfy the boundary conditions. This may involve boundary-condition treatment, interpolation, or other correction operations
    // reset counter of lost volume
    //if (verbose_) Pout << "Total particle volume neglected: " << tooMuch_<< endl; //Output the neglected void-fraction amount due to the setting above 0.1
    //tooMuch_ = 0.;
     if (verbose_)  	 
		 {
			 
			 
		 }
    // bring voidfraction from Eulerian Field to particle array
    //interpolationCellPoint<scalar> voidfractionInterpolator_(voidfractionNext_);
    //scalar voidfractionAtPos(0);

    for(int index=0; index< particleCloud_.numberOfParticles(); index++)
    {
        label cellID = particleCloud_.cellIDs()[index][0];

        if(cellID >= 0)
        {
            voidfractions[index][0] = voidfractionNext_[cellID];
        }
        else
        {
            voidfractions[index][0] = -1.;
        }
    }
  
}

inline double CGVoidFraction::particleRadius(label index) const
{
    return particleCloud_.radius(index);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

