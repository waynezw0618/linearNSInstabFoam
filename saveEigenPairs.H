
      volVectorField rEig
      (
          IOobject
          (
              "rEig"+name(i),
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionSet(0,1,-1,0,0,0,0),
          "zeroGradient"
      );

      volVectorField iEig
      (
          IOobject
          (
              "iEig"+name(i),
              runTime.timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::AUTO_WRITE
          ),
          mesh,
          dimensionSet(0,1,-1,0,0,0,0),
         "zeroGradient"
      );

      lts.copyEigenpairsToOpenFOAM(xr,rEig);
      lts.copyEigenpairsToOpenFOAM(xi,iEig);

      rEig.write();
      iEig.write();
