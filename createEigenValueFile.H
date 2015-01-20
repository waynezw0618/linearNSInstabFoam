    autoPtr<OFstream> EigenValueFilePtr_(NULL);

    if (Pstream::myProcNo() == 0)
    {
        if (EigenValueFilePtr_.empty())
        {
            word startTimeName =
                runTime.time().timeName(runTime.time().startTime().value());
            fileName EigenValueDir = runTime.time().path();

            // Create directory if does not exist.
            mkDir(EigenValueDir);

            // Open new file at start up
            EigenValueFilePtr_.reset(new OFstream(EigenValueDir/("EigenValue.dat")));

            // Add headers to output data
            if (EigenValueFilePtr_.valid())
            {
                EigenValueFilePtr_()
                    << "Solution method: " << type << nl
                    << "Number of requested eigenvalues: " << nev << nl
                    << "Stopping condition: tol=: " << tol << " maxit= " << maxit << nl
                    << "Stop at the " << its << "th iteration" << nl
                    << "Number of converged eigenpairs: "<< its << nl
                    << "ID " << "EigenValue " << "Error" << nl;
                EigenValueFilePtr_()<< endl;
            }
        }
    }
