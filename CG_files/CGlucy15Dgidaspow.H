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

    Gidaspow drag law
    - only valid for low-Reynolds number systems (Re_p<1000)
    - including interpolation of the velocity to the exact position
    - splits off explicit drag component due to fluctuation in fluid and particle
      velocity (optional via forceSubModel "ImExCorr")

Class
    CGlucy15Dgidaspow

SourceFiles
    CGlucy15Dgidaspow.C

\*---------------------------------------------------------------------------*/

#ifndef CGlucy15Dgidaspow_H
#define CGlucy15Dgidaspow_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
#include "voidFractionModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
const int numberOfMarkerPoints = 181;
/*---------------------------------------------------------------------------*\
                           Class CGlucy15Dgidaspow Declaration
\*---------------------------------------------------------------------------*/

class CGlucy15Dgidaspow
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;
	

    //mutable volScalarField voidfractioncgpp_;//add

	
	//const dimensionSet mydimension(0, 0, 0, 0, 0, 0, 0);//add
		
    const scalar phi_;

    word UsFieldName_;   

    const volVectorField& UsField_;  // the average particle velocity field

    mutable scalar switchingVoidfraction_; //voidfraction above which dilute formulation will be used
 
    //mutable scalar Voidfractioncgp_;// add by lyx
	
	vector offsets[numberOfMarkerPoints];// add by lyx
	
    virtual inline scalar Vp(int index, scalar radius) const
        {
            return 4.188790205*radius*radius*radius; //4/3*pi=4.188790205
        };
		
	//mutable double ** cellsPerParticle_;	
public:

    //- Runtime type information
    TypeName("CGlucy15Dgidaspow");


    // Constructors

        //- Construct from components
        CGlucy15Dgidaspow
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~CGlucy15Dgidaspow();


    // Member Functions
        void setForce() const;
		
	  //inline volScalarField& voidFractionm()const { return voidfractionm_; }//add by lyx
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
