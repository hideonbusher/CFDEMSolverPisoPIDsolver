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

#include "CGlucy15Dgidaspow.H"
#include "addToRunTimeSelectionTable.H"
//#include "averagingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(CGlucy15Dgidaspow, 0);

addToRunTimeSelectionTable
(
    forceModel,
    CGlucy15Dgidaspow,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
CGlucy15Dgidaspow::CGlucy15Dgidaspow
(
    const dictionary& dict,
    cfdemCloud& sm
)
:
    forceModel(dict,sm),
    propsDict_(dict.subDict(typeName + "Props")),
    velFieldName_(propsDict_.lookup("velFieldName")),
    U_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),
	//Urelate_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),//
	//Uparticlecgpp_(sm.mesh().lookupObject<volVectorField> (velFieldName_)),//
	//Jparticlecgpp_(sm.mesh(),dimensionedVector("zero", dimensionSet(1,1,-1,0,0), vector(0,0,0))), // 
	//roucgpp_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),//
	voidfractionFieldName_(propsDict_.lookup("voidfractionFieldName")),
    voidfraction_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),//
	//voidfractioncgpp_(sm.mesh().lookupObject<volScalarField> (voidfractionFieldName_)),//a
    phi_(propsDict_.lookupOrDefault<scalar>("phi",1.)),
    UsFieldName_(propsDict_.lookup("granVelFieldName")),
    UsField_(sm.mesh().lookupObject<volVectorField> (UsFieldName_)),
	switchingVoidfraction_(0.8)
{
    //Append the field names to be probed
    // suppress particle probe
    if (probeIt_ && propsDict_.found("suppressProbe"))
        probeIt_=!Switch(propsDict_.lookup("suppressProbe"));
    if(probeIt_)
    {
        particleCloud_.probeM().initialize(typeName, typeName+".logDat");
        particleCloud_.probeM().vectorFields_.append("dragForce"); //first entry must  be the force
        particleCloud_.probeM().vectorFields_.append("Urel");
        particleCloud_.probeM().scalarFields_.append("Rep");
        particleCloud_.probeM().scalarFields_.append("betaP");
        particleCloud_.probeM().scalarFields_.append("voidfraction");
        particleCloud_.probeM().writeHeader();
    }

    // init force sub model
    setForceSubModels(propsDict_);
    // define switches which can be read from dict
    forceSubM(0).setSwitchesList(0,true); // activate treatExplicit switch
    forceSubM(0).setSwitchesList(2,true); // activate implDEM switch
    forceSubM(0).setSwitchesList(3,true); // activate search for verbose switch
    forceSubM(0).setSwitchesList(4,true); // activate search for interpolate switch
    forceSubM(0).setSwitchesList(8,true); // activate scalarViscosity switch

    // read those switches defined above, if provided in dict
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).readSwitches();

    particleCloud_.checkCG(true);

    if (propsDict_.found("switchingVoidfraction"))
        switchingVoidfraction_ = readScalar(propsDict_.lookup("switchingVoidfraction"));  
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

CGlucy15Dgidaspow::~CGlucy15Dgidaspow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void CGlucy15Dgidaspow::setForce() const
{

    const volScalarField& rhoField = forceSubM(0).rhoField();
    #if defined(version24Dev)
       // there seems to have been a change in the return value of 
       // particleCloud_.turbulence().nu() used by forceSubM(0).nuField();
       const volScalarField& nufField = particleCloud_.turbulence().nu();
    #else
        const volScalarField& nufField = forceSubM(0).nuField();
    #endif


    //update force submodels to prepare for loop
    for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
        forceSubM(iFSub).preParticleLoop(forceSubM(iFSub).verbose());

    vector position(0,0,0);
    scalar voidfraction(1);
    vector Ufluid(0,0,0);
    vector drag(0,0,0);
    label cellI=0;

    vector Us(0,0,0);
    vector Ur(0,0,0);
    scalar ds(0);
    scalar dParcel(0);
    scalar nuf(0);
    scalar rho(0);
    scalar magUr(0);
    scalar Rep(0);
    scalar Vs(0);
    scalar localPhiP(0);
    scalar CdMagUrLag(0);       //Cd of the very particle
    scalar betaP(0);             //momentum exchange of the very particle
    scalar Cd(0);

    vector dragExplicit(0,0,0);
    scalar dragCoefficient(0);
	
	//roucgpp_ = 0.;             //
    //voidfractioncgpp_ = 1.;   //
	//Jparticlecgpp_ = Foam::vector(0, 0, 0);     //
	//Uparticlecgpp_ = Foam::vector(0, 0, 0);    //
	//Urelate_ = Foam::vector(0, 0, 0);         //
	

	scalar volume(0.);//add by lyx
	scalar radius(-1.);//add by lyx
	scalar voidfractionn(0);//origin new voidfraction
	vector Ufluidd(0,0,0);//origin new Ufluid
	vector Urr(0,0,0);//origin new Ufluid
	vector Upartt(0,0,0);
	vector calcufluid(0,0,0);//origin new Ufluid
	label i;//origin new index for cell

	label storesubcellkernal;
    scalar sumUweight; 
	scalar sumVweight; 
	
	scalar distance(0.);//new origin value
    scalar coefficient(0.);////new origin value
    scalar SumVol(0.);////new origin value
	scalar sumUCoefficient(0.);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    #include "resetVoidfractionInterpolator.H"
    #include "resetUInterpolator.H"
    #include "setupProbeModel.H"
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////voidfractioncgpp_////////////////////////////////////////////////////////////////////////////
//int np = particleCloud_.numberOfParticles();
//vector cellsperparticle_[np]; //define

//int storesubcell_[np][225];
//memset(storesubcell_, 0, sizeof(storesubcell_));
//int hang = 181;
//std::vector<std::vector<int>> storesubcell_(np, std::vector<int>(hang, 0.0));
//int text[np];
//int textt[np];

//int norepeat_[np][30];
//memset(norepeat_, 0, sizeof(norepeat_));





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(int index = 0;index <  particleCloud_.numberOfParticles(); ++index)
    {

			radius = particleCloud_.radius(index);//add by lyx 
			volume = Vp(index,radius);//add by lyx
			position = particleCloud_.position(index);
			Us = particleCloud_.velocity(index);//提取颗粒速度信息

			ds = 2*particleCloud_.radius(index);
            cellI = particleCloud_.cellIDs()[index][0];
            drag = vector(0,0,0);
            dragExplicit = vector(0,0,0);
            betaP = 0;
            Vs = 0;
			calcufluid =vector(0,0,0);
			voidfraction=0;
			voidfractionn=0;//add by lyx
			distance=0;//add by lyx
            dragCoefficient = 0;
            sumUCoefficient = 0; //////!!!!!!!!!!!!!!!!!!!!!!!!!!!which is lose in some other codes
			SumVol = 0; //////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!which is lose in some other codes
			coefficient= 0;
			Ufluidd =vector(0,0,0);
			Ufluid =vector(0,0,0);
			//Ufluid = U_[cellI];
			Upartt =vector(0,0,0);
			sumUweight=0;
			sumVweight=0;
            if (cellI > -1) // particle Found
            {
////////////////////////////////////////////////////////////////cal Ufluidd voidfractionn/////////////////////////////////////////////////////////
                Us = particleCloud_.velocity(index);
                ds = 2*particleCloud_.radius(index);
                dParcel = ds;
                forceSubM(0).scaleDia(ds,index); //caution: this fct will scale ds!
                rho = rhoField[cellI];
                nuf = nufField[cellI];		
				radius = particleCloud_.radius(index);//add by lyx 
				position = particleCloud_.position(index);
				voidfraction = voidfractionInterpolator_().interpolate(position,cellI);//debug			
////////////////////////////////////kernal1new////////////////////////////////////////////////////////////

			    for (int i = 0; i < particleCloud_.cellsPerParticle()[index][0]; i++) //kernalcell
				{
				storesubcellkernal = particleCloud_.cellIDs()[index][i];
				Ufluidd=Ufluidd+U_[storesubcellkernal]*particleCloud_.particleWeights()[index][i];
	            voidfractionn=voidfractionn+voidfraction_[storesubcellkernal]*particleCloud_.particleWeights()[index][i]*particleCloud_.mesh().V()[storesubcellkernal];
				sumUweight=sumUweight+particleCloud_.particleWeights()[index][i];
				sumVweight=sumVweight+particleCloud_.particleWeights()[index][i]*particleCloud_.mesh().V()[storesubcellkernal];
				}
				Ufluidd=Ufluidd/sumUweight;
				voidfractionn=voidfractionn/sumVweight;
				Urr=(Ufluidd)-(Us);

				
			
				if (voidfractionn < 0.05) voidfractionn = 0.05;
				if (voidfractionn > 1) voidfractionn = 1;
////////////////////////////////////kernal///////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                //Us = particleCloud_.velocity(index);
				//Ufluid = U_[cellI];
				//voidfraction=voidfraction_[cellI];
				

///////////////////Gidpsow///////////////////////////////////////////////////////////////////////////////////
                //Update any scalar or vector quantity
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                      forceSubM(iFSub).update(  index, 
                                                cellI,
                                                ds,
                                                Ufluid, 
                                                Us, 
                                                nuf,
                                                rho,
                                                forceSubM(0).verbose()
                                             );

                Vs = ds*ds*ds*M_PI/6;
                Ur = Urr;
				//Ur = Ufluidd-Us;
                magUr = mag(Ur);
                Rep=0.0;
                localPhiP = 1.0f-voidfractionn+SMALL;

                // calc particle's drag coefficient (i.e., Force per unit slip velocity and per m³ PARTICLE)
                if(voidfractionn > switchingVoidfraction_) //dilute
                {
                    Rep=ds*voidfractionn*magUr/(nuf+SMALL);
                    CdMagUrLag = (24.0*nuf/(ds*voidfractionn)) //1/magUr missing here, but compensated in expression for betaP!
                                 *(scalar(1.0)+0.15*Foam::pow(Rep, 0.687));

                    betaP = 0.75*(                                  //this is betaP = beta / localPhiP!
                                            rho*voidfractionn*CdMagUrLag
                                          /
                                            (ds*Foam::pow(voidfractionn,2.65))
                                          );
                }
                else  //dense
                {
                    betaP = (150 * localPhiP*nuf*rho)          //this is betaP = beta / localPhiP!
                             /  (voidfractionn*ds*phi_*ds*phi_)
                            +
                              (1.75 * magUr * rho)
                             /((ds*phi_));
                }

                // calc particle's drag
                dragCoefficient = Vs*betaP;
                if (modelType_=="B")
                    dragCoefficient /= voidfractionn;

                forceSubM(0).scaleCoeff(dragCoefficient,dParcel,index);
                drag = dragCoefficient * Ur;

                // explicitCorr
                for (int iFSub=0;iFSub<nrForceSubModels();iFSub++)
                    forceSubM(iFSub).explicitCorr( drag, 
                                                   dragExplicit,
                                                   dragCoefficient,
                                                   Ufluid, U_[cellI], Us, UsField_[cellI],
                                                   forceSubM(iFSub).verbose()
                                                 );
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                   /* if(forceSubM(0).verbose() && (index >=6000 && index <6003 || index >=700 && index <703))
                {
					//Pout << "cellI = " << cellI << endl;
                    //Pout << "index = " << index << endl;
					//Pout << "text[index] = " << text[index] << endl;
                }  */

                //Set value fields and write the probe
                if(probeIt_)
                {
                    #include "setupProbeModelfields.H"
                    // Note: for other than ext one could use vValues.append(x)
                    // instead of setSize
                    vValues.setSize(vValues.size()+1, drag);           //first entry must the be the force
                    vValues.setSize(vValues.size()+1, Ur);
                    sValues.setSize(sValues.size()+1, Rep); 
                    sValues.setSize(sValues.size()+1, betaP);
                    sValues.setSize(sValues.size()+1, voidfraction);
                    particleCloud_.probeM().writeProbe(index, sValues, vValues);
                }
            }

            // write particle based data to global array
            forceSubM(0).partToArray(index,drag,dragExplicit,Ufluid,dragCoefficient);

        //}// end if mask
    }// end loop particles
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
