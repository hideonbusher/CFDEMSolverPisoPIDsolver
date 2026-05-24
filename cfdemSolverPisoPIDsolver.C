/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright (C) 1991-2009 OpenCFD Ltd.
                                Copyright (C) 2009-2012 JKU, Linz
                                Copyright (C) 2012-     DCS Computing GmbH,Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling.  If not, see <http://www.gnu.org/licenses/>.

Application
    cfdemSolverPiso

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.
    The code is an evolution of the solver pisoFoam in OpenFOAM(R) 1.6,
    where additional functionality for CFD-DEM coupling is added.
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"

#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif
#if defined(versionv1606plus) || defined(version40)
    #include "fvOptions.H"
#else
    #include "fvIOoptionList.H"
#endif
#include "fixedFluxPressureFvPatchScalarField.H"
#include "cfdemCloud.H"

#if defined(anisotropicRotation)
    #include "cfdemCloudRotation.H"
#endif
#if defined(superquadrics_flag)
    #include "cfdemCloudRotationSuperquadric.H"
#endif
#include "implicitCouple.H"
#include "clockModel.H"
#include "smoothingModel.H"
#include "forceModel.H"

#include <iomanip>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #if defined(version30)
        pisoControl piso(mesh);
        #include "createTimeControls.H"
    #endif
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    // create cfdemCloud
    #include "readGravitationalAcceleration.H"
    #include "checkImCoupleM.H"
    #if defined(anisotropicRotation)
        cfdemCloudRotation particleCloud(mesh);
    #elif defined(superquadrics_flag)
        cfdemCloudRotationSuperquadric particleCloud(mesh);
    #else
        cfdemCloud particleCloud(mesh);
    #endif
    #include "checkModelType.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #if defined(version30)
            #include "readTimeControls.H"
            #include "CourantNo.H"
            #include "setDeltaT.H"
        #else
            #include "readPISOControls.H"
            #include "CourantNo.H"
        #endif

        // do particle stuff
        particleCloud.clockM().start(1,"Global");
        particleCloud.clockM().start(2,"Coupling");
        bool hasEvolved = particleCloud.evolve(voidfraction,Us,U);

        if(hasEvolved)
        {
            particleCloud.smoothingM().smoothenAbsolutField(particleCloud.forceM(0).impParticleForces());
        }
    
        Ksl = particleCloud.momCoupleM(particleCloud.registryM().getProperty("implicitCouple_index")).impMomSource();
        Ksl.correctBoundaryConditions();

        surfaceScalarField voidfractionf = fvc::interpolate(voidfraction);
        phi = voidfractionf*phiByVoidfraction;

        //Force Checks
        #include "forceCheckIm.H"

        #include "solverDebugInfo.H"
        particleCloud.clockM().stop("Coupling");

        particleCloud.clockM().start(26,"Flow");
        
static scalar error_integral = 0.0;
static scalar Q_error_last = 0.0;  // Used for D control
	
// Get the patch ID of the right boundary
//label rightPatchID = mesh.boundaryMesh().findPatchID("right");
// Calculate the current flow rate Q_current at the right boundary
//scalar Q_current = gSum(mag(phi.boundaryField()[rightPatchID])); 
//scalar Q_current = gSum((phi.boundaryField()[rightPatchID])); 

scalar Q_current = 0.0;

// Custom target profile location x = x_plane
//scalar x_plane = 0.1;

// Small tolerance to avoid floating-point comparison issues
scalar tol = 1e-6;

const auto& phiField = phi.internalField();  
const vectorField& Cf = mesh.Cf();  // Face-center coordinates of all faces // Cf is a vectorField representing the center coordinate of each face
const vectorField& Sf = mesh.Sf();  // Area vectors of all faces (direction + magnitude)

forAll(phiField, faceI)
{
    // Select internal faces only
    if (faceI < mesh.nInternalFaces())
    {
        if (mag(Cf[faceI].x() - x_plane) < tol)
        {
            
            Q_current += phiField[faceI];
        }
    }
}

// Parallel reduction summation
reduce(Q_current, sumOp<scalar>());

// PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
scalar Q_error = Q_target - Q_current;


// IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
// scalar max_integral = 10.0; // Set the maximum integral value to avoid excessive accumulation
// **Add an anti-windup mechanism to prevent an excessively large integral term** relative error tolerance
scalar error_ratio = mag(Q_error) / (mag(Q_target) + SMALL);
if (error_ratio < 0.0025)
{
    error_integral *= 0.9999;  // Gradually decay the integral term
}
else
{
    error_integral += Q_error * min(runTime.deltaT().value(), 0.01);
}
// Limit the integral error to prevent over-compensation
error_integral = min(max(error_integral, -max_integral), max_integral);


//DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
scalar dt = runTime.deltaT().value();
scalar Q_error_derivative = (Q_error - Q_error_last) / (dt + SMALL);
Q_error_last = Q_error;  

// Calculate the P, I, and D control terms
scalar P_term = Kp * Q_error;
scalar I_term = Ki * error_integral;
scalar D_term = Kd * Q_error_derivative;
// Calculate the volume-force correction term (including I control)
scalar Fvol_scalar = P_term + I_term+D_term;

if (Pstream::master() && runTime.timeIndex() % 10 == 0)
{
    std::ofstream pidLogFile;
    if (runTime.timeIndex() == 0)
    {
        pidLogFile.open("PID_log.csv");
        pidLogFile << "Time,Q_current,Q_error,error_integral,P,I,D,Fvol_scalar\n";
    }
    else
    {
        pidLogFile.open("PID_log.csv", std::ios_base::app);
    }

    pidLogFile << std::fixed << std::setprecision(6)
               << runTime.timeName() << ","
               << Q_current << ","
               << Q_error << ","
               << error_integral << ","
               << P_term << ","
               << I_term << ","
               << D_term << ","
               << Fvol_scalar << "\n";

    pidLogFile.close();
}


// Volume force applied in the x direction
dimensionedVector Fvol
(
    "Fvol",   // Variable name
    dimensionSet(0,1,-2,0,0,0,0), // Unit of volume force [m/s²]
    vector(Fvol_scalar, 0, 0)     // Force in the x direction
);


		
        if(particleCloud.solveFlow())
        {
            // Pressure-velocity PISO corrector
            {
                // Momentum predictor
                fvVectorMatrix UEqn
                (
                    fvm::ddt(voidfraction,U) - fvm::Sp(fvc::ddt(voidfraction),U)
                  + fvm::div(phi,U) - fvm::Sp(fvc::div(phi),U)
//                + turbulence->divDevReff(U)
                  + particleCloud.divVoidfractionTau(U, voidfraction)
                  ==
                  - fvm::Sp(Ksl/rho,U)
                  + fvOptions(U)
				  + gravity*voidfraction  //add by Liuyx
				  + drive*voidfraction  //add by Liuyx
				  + Fvol*voidfraction
                );

                UEqn.relax();
                fvOptions.constrain(UEqn);

                #if defined(version30)
                    if (piso.momentumPredictor())
                #else
                    if (momentumPredictor)
                #endif
                {
                    if (modelType=="B" || modelType=="Bfull")
                        solve(UEqn == - fvc::grad(p) + Ksl/rho*Us);
                    else
                        solve(UEqn == - voidfraction*fvc::grad(p) + Ksl/rho*Us);

                    fvOptions.correct(U);
                }

                // --- PISO loop   PISO pressure Poisson equation solution
                #if defined(version30)
                    while (piso.correct())
                #else
                    for (int corr=0; corr<nCorr; corr++)
                #endif
                {
                    volScalarField rUA = 1.0/UEqn.A();

                    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));
                    volScalarField rUAvoidfraction("(voidfraction2|A(U))",rUA*voidfraction);
                    surfaceScalarField rUAfvoidfraction("(voidfraction2|A(U)F)", fvc::interpolate(rUAvoidfraction));

                    U = rUA*UEqn.H();

                    #ifdef version23
                        phi = ( fvc::interpolate(U) & mesh.Sf() )
                            + rUAfvoidfraction*fvc::ddtCorr(U, phiByVoidfraction);
                    #else
                        phi = ( fvc::interpolate(U) & mesh.Sf() )
                            + fvc::ddtPhiCorr(rUAvoidfraction, U, phiByVoidfraction);
                    #endif
                    surfaceScalarField phiS(fvc::interpolate(Us) & mesh.Sf());
                    phi += rUAf*(fvc::interpolate(Ksl/rho) * phiS);

                    if (modelType=="A")
                        rUAvoidfraction = volScalarField("(voidfraction2|A(U))",rUA*voidfraction*voidfraction);

                    // Update the fixedFluxPressure BCs to ensure flux consistency
                    #include "fixedFluxPressureHandling.H"
                    

                    // Non-orthogonal pressure corrector loop for non-orthogonal mesh correction
                    #if defined(version30)
                        while (piso.correctNonOrthogonal())
                    #else
                        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                    #endif
                    {
                        // Pressure corrector: solution of the pressure Poisson equation
                        fvScalarMatrix pEqn
                        (
                            fvm::laplacian(rUAvoidfraction, p) == fvc::div(voidfractionf*phi) + particleCloud.ddtVoidfraction()
                        );
                        pEqn.setReference(pRefCell, pRefValue);

                        #if defined(version30)
                            pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
                            if (piso.finalNonOrthogonalIter())    // Non-orthogonal mesh pressure correction
                            {
                                phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                            }
                        #else
                            if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                                #if defined(versionExt32)
                                    pEqn.solve(mesh.solutionDict().solver("pFinal"));
                                #else
                                    pEqn.solve(mesh.solver("pFinal"));
                                #endif
                            else
                                pEqn.solve();

                            if (nonOrth == nNonOrthCorr)
                            {
                                phiByVoidfraction = phi - pEqn.flux()/voidfractionf;
                            }
                        #endif

                    } // end non-orthogonal corrector loop for non-orthogonal mesh correction

                    phi = voidfractionf*phiByVoidfraction; // Correct the face flux
                    #include "continuityErrorPhiPU.H"

                    if (modelType=="B" || modelType=="Bfull")
                        U -= rUA*fvc::grad(p) - Ksl/rho*Us*rUA;
                    else
                        U -= voidfraction*rUA*fvc::grad(p) - Ksl/rho*Us*rUA;

                    U.correctBoundaryConditions();  // Correct the velocity and velocity boundary conditions
                    fvOptions.correct(U);

                } // end piso loop
            }

            laminarTransport.correct();
            turbulence->correct();
        }// end solveFlow
        else
        {
            Info << "skipping flow solution." << endl;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        particleCloud.clockM().stop("Flow");
        particleCloud.clockM().stop("Global");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

