/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::linearizedTimeStepper

Description
    Linear stability analysis

SourceFiles
    linearizedTimeStepper.C

\*---------------------------------------------------------------------------*/

#ifndef linearizedTimeStepper_H
#define linearizedTimeStepper_H

#include <slepceps.h>
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Matrix wrapper
PetscErrorCode MatVecMult(Mat A, Vec x, Vec y);

/*---------------------------------------------------------------------------*\
                           Class linearizedTimeStepper Declaration
\*---------------------------------------------------------------------------*/

class linearizedTimeStepper
{
    // Private data
        Time *runTimePtr;
        fvMesh *meshPtr;
        volScalarField& p;
        volVectorField& U;
        volVectorField& Up;
        surfaceScalarField& phi;
        label pRefCell;
        label loopID;
        label cellNumber;
        label N;
        scalar pRefValue;
        scalar cumulativeContErr;
        dimensionedScalar endTimeValue;
        word timeStepper;

    // Private Member Functions
        // Copy SLEPc data to OpenFOAM field
        void copySLEPcDataToOpenFOAMField(Vec x);

        // Copy OpenFOAM Field to SLEPc data
        void copyOpenFOAMFieldToSLEPcData(Vec y);

        // PISO method
        void pisoStepper();

        // Chorin method
        void chorinStepper();

        // RK4 method
        void rk4Stepper();

public:

    //- Runtime type information
    TypeName("linearizedTimeStepper");

    // Constructors

        //- Construct from components
        linearizedTimeStepper
        (
            Time   &runTime_,
            fvMesh &mesh_,
            label pRefCell_,
            scalar pRefValue_,
            word timeStepper_
        );

    //- Destructor
    virtual ~linearizedTimeStepper()
    {}


    // Member Functions

        // Matrix free code
        PetscErrorCode MatVecMult(Vec x, Vec y);

        // Copy final eigenpairs to OpenFOAM and save the data
        void copyEigenpairsToOpenFOAM(Vec v, volVectorField &ev);

        // Set initial guess
        void setInitialVector(Vec y);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
