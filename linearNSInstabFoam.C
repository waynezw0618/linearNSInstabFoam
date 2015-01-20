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

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "linearizedTimeStepper.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initSlepc.H"

    //  Compute the operator matrix that defines the eigensystem, Ax=kx
    ierr = MatCreateShell(PETSC_COMM_WORLD,n,n,PETSC_DETERMINE,PETSC_DETERMINE,NULL,&A);CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);CHKERRQ(ierr);

    ierr = MatShellSetContext(A, &lts); CHKERRQ(ierr);
    ierr = MatShellSetOperation(A, MATOP_MULT, (void(*)()) MatVecMult); CHKERRQ(ierr);

    ierr = MatGetVecs(A, PETSC_NULL,&xr); CHKERRQ(ierr);
    ierr = MatGetVecs(A, PETSC_NULL,&xi); CHKERRQ(ierr);
    ierr = MatGetVecs(A, PETSC_NULL,&iv); CHKERRQ(ierr);

    //  Create the eigensolver and set various options
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

    //  Set operators. In this case, it is a standard eigenvalue problem
    ierr = EPSSetOperators(eps, A, PETSC_NULL);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_NHEP);CHKERRQ(ierr);
    ierr = EPSSetDimensions(eps, nev, ncv, mpd);CHKERRQ(ierr);
    ierr = EPSSetWhichEigenpairs(eps,which);CHKERRQ(ierr);
    ierr = EPSSetType(eps,method.c_str());CHKERRQ(ierr);
    ierr = EPSSetTolerances(eps,s_tol,s_maxit);CHKERRQ(ierr);
    ierr = EPSSetConvergenceTest(eps,cov);CHKERRQ(ierr);

    //  Set solver parameters at runtime
    ierr = EPSSetFromOptions(eps); CHKERRQ(ierr);

    //  Set initial vector
    lts.setInitialVector(iv);
    ierr = EPSSetInitialSpace(eps,1,&iv);CHKERRQ(ierr);

    ierr = EPSView(eps, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    //  Solve the eigensystem
    ierr = EPSSolve(eps);CHKERRQ(ierr);

    //  Optional: Get some information from the solver and display it
    ierr = EPSGetType(eps,&type);CHKERRQ(ierr);

    ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4G, maxit=%D\n",tol,maxit);CHKERRQ(ierr);
    ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD," Stop at the %Dth iterations\n",its);CHKERRQ(ierr);

    //  Get number of converged approximate eigenpairs
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
    // ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);
    #include "createEigenValueFile.H"
    if (nconv>0)
    {
        // Display eigenvalues and relative errors
        for (i=0;i<nev;i++)
        {
            /*
             *         Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
             *                 ki (imaginary part)
             */
             ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
             /*
              *          Compute the relative error associated to each eigenpair
              */
             //ierr = EPSComputeRelativeError(eps,i,&error);CHKERRQ(ierr);

             #if defined(PETSC_USE_COMPLEX)
                 re = PetscRealPart(kr);
                 im = PetscImaginaryPart(kr);
             #else
                 re = kr;
                 im = ki;
             #endif

             #include "saveEigenPairs.H"

             if (Pstream::myProcNo() == 0)
             {
                 if (im!=0.0)
                 {
                     //ierr = PetscPrintf(PETSC_COMM_WORLD," %9F%+9F j %12G\n",re,im,error);CHKERRQ(ierr);
                     EigenValueFilePtr_()<< i << tab << re <<"+"<< im << "j"<<tab << error << nl;
                 }
                 else
                 {
                     //ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12F       %12G\n",re,error);CHKERRQ(ierr);
                     EigenValueFilePtr_()<< i << tab << re << tab << error << nl;
                 }
                 EigenValueFilePtr_()<< endl;
             }
        }
        //ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
    }

    //ierr = EPSPrintSolution(eps,PETSC_NULL);CHKERRQ(ierr);
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
    ierr = MatDestroy(&A);CHKERRQ(ierr);
    ierr = SlepcFinalize();CHKERRQ(ierr);
    return (0);
}

// ************************************************************************* //
