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

Class
    CGVoidFraction

SourceFiles
    CGVoidFraction.C

\*---------------------------------------------------------------------------*/

#ifndef CGVoidFraction_H
#define CGVoidFraction_H
#include "voidFractionModel.H"
#include "vectorField.H" // 用于处理 vectorField 类型
#include "word.H"
#include "fvMesh.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
const int numberOfMarkerPoints = 313;

/*---------------------------------------------------------------------------*\
                           Class noDrag Declaration
\*---------------------------------------------------------------------------*/

class CGVoidFraction
:
    public voidFractionModel
{

private:
        dictionary propsDict_;

        const fvMesh& mesh_;

       // label nX_, nY_, nZ_; // 
   
       // mutable volVectorField A;
   
        //IOdictionary blockMeshDict_;
   
        bool verbose_;

        Switch procBoundaryCorrection_;
        Switch periodicx_;
		Switch periodicy_;
		Switch periodicz_;
		Switch TopologyAccel_;
        const scalar alphaMin_;          // min value of voidFraction

	    mutable bool alphaLimited_;

        mutable scalar tooMuch_;         // particle volume which is lost due to voidFraction limitation

        bool interpolation_;

        bool cfdemUseOnly_;

        virtual inline scalar Vp(int index, scalar radius, scalar scaleVol) const
        {
            return 4.188790205*radius*radius*radius*scaleVol; //4/3*pi=4.188790205
        };

        vector offsets[numberOfMarkerPoints];

        const scalar meshx_; const scalar meshy_; const scalar meshz_; const scalar decx_; const scalar decy_; const scalar decz_;
		label calculateLocalGridID(
         const vector& subPosition,
         const Foam::boundBox& globalBb,
         int meshx, int meshy, int meshz,
         Foam::scalar decx_, Foam::scalar decy_, Foam::scalar decz_) const;
public:

    //- Runtime type information
    TypeName("CG");


    // Constructors

        //- Construct from components
        CGVoidFraction
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~CGVoidFraction();


    // Member Functions
        void setvoidFraction(double** const& ,double**&, double**&, double**&, double**&) const;

        virtual double particleRadius(label) const;
		
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
