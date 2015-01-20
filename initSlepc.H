    Mat            A;               /* operator matrix */
    EPS            eps;             /* eigenproblem solver context */
    const EPSType  type;
    PetscMPIInt    size;
    PetscInt       i,n,mpd=100,nconv,its,maxit,s_maxit,nev,ncv;
    PetscReal      error=0,tol,re, im, s_tol;
    PetscScalar    kr,ki;
    Vec            xr,xi,iv;
    PetscErrorCode ierr;
    EPSConv        cov;
    word           method;
    EPSWhich       which;
    word           TimeStepper;

    {
        IOdictionary slepcParameters
        (
            IOobject
            (
                "slepcParameters",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        List<string> slepcArgs(slepcParameters.lookup("arguments"));
        method=slepcParameters.lookupOrDefault<word>("EigenProblemSolver","krylovschur");
        which=(EPSWhich) slepcParameters.lookupOrDefault<PetscInt>("ProblemType",1);
        cov=(EPSConv) slepcParameters.lookupOrDefault<PetscInt>("ConvergenceCriterion",1);
        s_tol=slepcParameters.lookupOrDefault<PetscReal>("Tolerance",1.0e-6);
        ncv=slepcParameters.lookupOrDefault<PetscInt>("KrylovspaceSize",8);
        nev=slepcParameters.lookupOrDefault<PetscInt>("NumEigenVaules",4);
        s_maxit=slepcParameters.lookupOrDefault("MaxIterationNumber",10000);
        TimeStepper=slepcParameters.lookupOrDefault("TimeStepper",word("piso"));

        Info<< "The Torelence is set to "<<s_tol<<nl
            << "The Max IterationNumber is set to "<<s_maxit<<nl
            << "The Krylov space size is set to "<<ncv<<nl
            << "The Number of requested EigenVaules is set to "<<nev<<nl
            <<endl;

        int slepcArgc = slepcArgs.size() + 1;
        char **slepcArgv = NULL;

        char progName [] = "LineariseNSFoam";
        slepcArgv = new char* [slepcArgc];
        int PETSC_MAX_LINE_LENGTH = 32768;
        slepcArgv[0] = new char [PETSC_MAX_LINE_LENGTH];
        sprintf(slepcArgv[0], "%s", progName);
        forAll(slepcArgs, argI)
        {
            slepcArgv[argI + 1] = new char [ PETSC_MAX_LINE_LENGTH ];
            sprintf(slepcArgv[argI+1], "%s", slepcArgs[argI].c_str());
        }
        SlepcInitialize(&slepcArgc, &slepcArgv, PETSC_NULL, PETSC_NULL);

        delete [] slepcArgv[0];
        forAll(slepcArgs, argI)
        {
            delete [] slepcArgv[argI + 1];
        }
        delete [] slepcArgv;

    }

    linearizedTimeStepper lts(runTime, mesh, pRefCell, pRefValue,TimeStepper);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    n=2*mesh.nCells();
