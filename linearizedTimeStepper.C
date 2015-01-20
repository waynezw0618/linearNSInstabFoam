/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "linearizedTimeStepper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(linearizedTimeStepper, 0);


PetscErrorCode MatVecMult(Mat A, Vec x, Vec y)
{
    linearizedTimeStepper *lts;
    MatShellGetContext(A, &lts);
    return lts->MatVecMult(x, y);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void linearizedTimeStepper::copySLEPcDataToOpenFOAMField(Vec x)
{
    PetscScalar *px;
    VecGetArray(x,&px);

    fvMesh &mesh = *meshPtr;
    forAll(mesh.C(),celli)
    {
        Up.internalField()[celli].x()=px[celli];
        Up.internalField()[celli].y()=px[celli+cellNumber];
        //Up.internalField()[celli].z()=x[celli+2*cellNumber];
    }
};


void linearizedTimeStepper::copyOpenFOAMFieldToSLEPcData(Vec y)
{
    PetscInt istart;
    VecGetOwnershipRange(y,&istart,PETSC_NULL);

    fvMesh &mesh = *meshPtr;
    forAll(mesh.C(),celli)
    {
        int i = celli + istart;
        int j = i + cellNumber;
        VecSetValues(y,1,&i,&Up.internalField()[celli].x(),INSERT_VALUES);
        VecSetValues(y,1,&j,&Up.internalField()[celli].y(),INSERT_VALUES);
    }

    VecAssemblyBegin(y);
    VecAssemblyEnd(y);
};

void linearizedTimeStepper::pisoStepper()
{
    Info<< "\nStarting "<< loopID <<"th time Stepper loop\n"   << endl;

    fvMesh &mesh = *meshPtr;
    Time &runTime= *runTimePtr;

    incompressible::turbulenceModel &turbulence
        = const_cast<incompressible::turbulenceModel&>(mesh.lookupObject<incompressible::turbulenceModel>("turbulenceModel"));

    cumulativeContErr = 0;

    #include "createUpphi.H"

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");

        const int nCorr =
            pisoDict.lookupOrDefault<int>("nCorrectors", 1);

        const int nNonOrthCorr =
            pisoDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

        const bool momentumPredictor =
            pisoDict.lookupOrDefault("momentumPredictor", true);

        #include "CourantNo.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(Up)
              + fvm::div(phi, Up)
              + turbulence.divDevReff(Up)
            );

            UEqn.relax();

            if (momentumPredictor)
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop

            for (int corr=0; corr<nCorr; corr++)
            {
                volScalarField rAU(1.0/UEqn.A());

                Up = rAU*UEqn.H();

                upphi = (fvc::interpolate(Up) & mesh.Sf())
                      + fvc::ddtPhiCorr(rAU, Up, upphi);

                adjustPhi(upphi, Up, p);

                // Non-orthogonal pressure corrector loop
                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(upphi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        upphi -= pEqn.flux();
                    }
                }

                {
                    volScalarField contErr(fvc::div(upphi));

                    scalar sumLocalContErr = runTime.deltaTValue()*
                        mag(contErr)().weightedAverage(mesh.V()).value();

                    scalar globalContErr = runTime.deltaTValue()*
                        contErr.weightedAverage(mesh.V()).value();
                    cumulativeContErr += globalContErr;

                    Info<< "time step continuity errors : sum local = " << sumLocalContErr
                        << ", global = " << globalContErr
                        << ", cumulative = " << cumulativeContErr
                        << endl;
                }

                Up -= rAU*fvc::grad(p);
                Up.correctBoundaryConditions();
            }
        }
        turbulence.correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    }

    runTime.setTime(runTime.timeOutputValue(), runTime.timeIndex());
    runTime.setEndTime(runTime.endTime()+endTimeValue);

    Info<< "\nEnd of "<< loopID <<"th time stepper loop\n"   << endl;

    loopID++;
}

void linearizedTimeStepper::chorinStepper()
{
    Info<< "\nStarting "<< loopID <<"th time Stepper loop\n"   << endl;
    fvMesh &mesh = *meshPtr;
    Time &runTime= *runTimePtr;

    #include "createUpphi.H"

    incompressible::turbulenceModel &turbulence
        = const_cast<incompressible::turbulenceModel&>(mesh.lookupObject<incompressible::turbulenceModel>("turbulenceModel"));

    cumulativeContErr = 0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(Up)
          + fvm::div(phi, Up)
          - fvm::laplacian(turbulence.nuEff(), Up)
        );

        UEqn.solve();

        upphi = (fvc::interpolate(Up) & mesh.Sf());

        fvScalarMatrix pEqn
        (
            fvm::laplacian(p) == 1/runTime.deltaT()*fvc::div(upphi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        Up -= runTime.deltaT()*fvc::grad(p);
        Up.correctBoundaryConditions();


        {
            volScalarField contErr(fvc::div(upphi));

            scalar sumLocalContErr = runTime.deltaTValue()*
                mag(contErr)().weightedAverage(mesh.V()).value();

            scalar globalContErr = runTime.deltaTValue()*
                contErr.weightedAverage(mesh.V()).value();
            cumulativeContErr += globalContErr;

            Info<< "time step continuity errors : sum local = " << sumLocalContErr
                << ", global = " << globalContErr
                << ", cumulative = " << cumulativeContErr
                << endl;
        }

        turbulence.correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    runTime.setTime(runTime.timeOutputValue(), runTime.timeIndex());
    runTime.setEndTime(runTime.endTime()+endTimeValue);

    Info<< "\nEnd of "<< loopID <<"th time stepper loop\n"   << endl;

    loopID++;
}

void linearizedTimeStepper::rk4Stepper()
{
    Info<< "\nStarting "<< loopID <<"th time Stepper loop\n"   << endl;
    fvMesh &mesh = *meshPtr;
    Time &runTime= *runTimePtr;

    #include "createUpphi.H"

    incompressible::turbulenceModel &turbulence
        = const_cast<incompressible::turbulenceModel&>(mesh.lookupObject<incompressible::turbulenceModel>("turbulenceModel"));

    scalar a[4] = {1./6., 1./3., 1./3., 1./6.}, b[4] = {0.5, 0.5, 1, 1};

    cumulativeContErr = 0;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        dimensionedScalar dt("dt",dimensionSet(0, 0, 1, 0, 0), runTime.deltaTValue());

        #include "CourantNo.H"
        volVectorField Upold(Up),Uptemp(Up);

        for (int i=0;i<4; i++)
        {
            if (i > 0)
            {
                upphi = (fvc::interpolate(Up) & mesh.Sf());
            }
            volVectorField dUp = dt*(-fvc::div(phi, Up) + turbulence.nuEff()*fvc::laplacian(Up));
            if(i<3){Up = Upold + b[i]*dUp;}
            U.correctBoundaryConditions();

            fvScalarMatrix pEqn
            (
                fvm::laplacian(p) == fvc::div(Up)/dt
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            Up -= dt*fvc::grad(p);Up.correctBoundaryConditions();

            Uptemp+=a[i]*dUp;

            {
                volScalarField contErr(fvc::div(upphi));

                scalar sumLocalContErr = runTime.deltaTValue()*
                    mag(contErr)().weightedAverage(mesh.V()).value();

                scalar globalContErr = runTime.deltaTValue()*
                    contErr.weightedAverage(mesh.V()).value();
                cumulativeContErr += globalContErr;

                Info<< "time step continuity errors : sum local = " << sumLocalContErr
                    << ", global = " << globalContErr
                    << ", cumulative = " << cumulativeContErr
                    << endl;
            }
        }

        Up=Uptemp;
        Up = Up -dt*fvc::grad(p);Up.correctBoundaryConditions();

        turbulence.correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    runTime.setTime(runTime.timeOutputValue(), runTime.timeIndex());
    runTime.setEndTime(runTime.endTime()+endTimeValue);

    Info<< "\nEnd of "<< loopID <<"th time stepper loop\n"   << endl;

    loopID++;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

linearizedTimeStepper::linearizedTimeStepper
(
    Time   &runTime_,
    fvMesh &mesh_,
    label pRefCell_,
    scalar pRefValue_,
    word timeStepper_
)
:
    runTimePtr(&runTime_),
    meshPtr(&mesh_),
    p(const_cast<volScalarField&>(mesh_.lookupObject<volScalarField>("p"))),
    U(const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("U"))),
    Up(const_cast<volVectorField&>(mesh_.lookupObject<volVectorField>("Up"))),
    phi(const_cast<surfaceScalarField&>(mesh_.lookupObject<surfaceScalarField>("phi"))),
    pRefCell(pRefCell_),
    loopID(0),
    cellNumber(mesh_.nCells()),
    N(2*cellNumber),
    pRefValue(pRefValue_),
    cumulativeContErr(0),
    endTimeValue((*runTimePtr).endTime()),
    timeStepper(timeStepper_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

PetscErrorCode linearizedTimeStepper::MatVecMult(Vec x, Vec y)
{
    PetscFunctionBegin;

    copySLEPcDataToOpenFOAMField(x);

    if (timeStepper.compare("piso") == 0)
    {
        Info<< "The time stepping method is set to " << timeStepper << nl << endl;

        this->pisoStepper();
    }
    else if (timeStepper.compare("rk4") == 0)
    {
        Info<< "The time stepping method is set to " << timeStepper << nl << endl;
        this->rk4Stepper();
    }
    else if (timeStepper.compare("chorin") == 0)
    {
        Info<< "The time stepping method is set to " << timeStepper << nl << endl;
        this->chorinStepper();
    }
    else
    {
        Info<< "Please check that the time stepping method is correctly specified" << nl << endl;
    }

    copyOpenFOAMFieldToSLEPcData(y);

    PetscFunctionReturn(0);
}

void linearizedTimeStepper::copyEigenpairsToOpenFOAM(Vec v, volVectorField &ev)
{
    PetscScalar *pv;
    VecGetArray(v,&pv);

    fvMesh &mesh = *meshPtr;
    forAll(mesh.C(), celli)
    {
        ev.internalField()[celli].x() = pv[celli];
        ev.internalField()[celli].y() = pv[celli+cellNumber];
        ev.internalField()[celli].z() = 0;//pv[celli+2*cellNumber];
    }

    ev.correctBoundaryConditions();
};

void linearizedTimeStepper::setInitialVector(Vec y)
{
    copyOpenFOAMFieldToSLEPcData(y);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
