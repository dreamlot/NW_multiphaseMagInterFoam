    MRF.correctBoundaryVelocity(U);

    //NW 2017-11-13 21:11
    Info<< "Calculating solveEq and NOTsolveEq...\n" << endl;
    volScalarField  solveEq(pos(mixture.isSolid()-0.8));
    volScalarField  NOTsolveEq(1-solveEq);
    //volScalarField  isSolid(mixture.isSolid());
    //isSolid.write();

    //NW 2017-11-13 16:11
    volScalarField  tmpCoeff(
        IOobject
        (
            "tmpCoeff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tmpCoeff", dimDensity/dimTime, 1.0)
    );


/*
    //NW 2017-11-13 21:11
    volScalarField tmpCoeff("tmpCoeff", dimensionedScalar("tmpCoeff", dimDensity/dimTime, 1.0) );
*/
/*
    fvVectorMatrix UEqn
    (//NW 2017-11-12 14:11
        solveEq*
        (
            fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
          + MRF.DDt(rho, U)
          + turbulence->divDevRhoReff(rho, U)

//         ==
//            fvOptions(rho, U)

        )
        +
        NOTsolveEq*tmpCoeff*
        (
            fvm::Sp(0.0,U)
        )

    );
*/
    //NW 2017-11-13 21:11
    Info<< "Setting UEqn...\n" << endl;

    fvVectorMatrix UEqn
    (
        //NW 2017-11-16 12:35
        solveEq * (
            (
            fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
          + MRF.DDt(rho, U)
          + turbulence->divDevRhoReff(rho, U)
            )
         //==
         //NW 2017-11-16 12:36
         -
         (
           fvOptions(rho, U)
          +
            fvc::reconstruct
            (
                (
                   mixture.surfaceTensionForce()
                 - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh + p_m)
                  + constant::electromagnetic::mu0
                  * fvc::interpolate(mag(H))*fvc::snGrad(mag(H))
                ) * mesh.magSf()
            )
        )
        )
        //NW 2017-11-16 12:35
        +
        NOTsolveEq*tmpCoeff*
        (
            fvm::Sp(0.0,U)
        )
        //NW 2017-11-16 13:00
        //==0
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    /*
    //NW20170926
    Foam::tmp<Foam::surfaceScalarField> magneticNormalForce() const
    {
        tmp<volVectorField> tmnf
        (
            new volVectorField
            (
                IOobject
                (
                    "magneticNormalForce",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedVector
                (
                    "magneticNormalForce",
                    dimensionSet(1, -2, -2, 0, 0,0,0),
                    vector(0.0,0.0,0.0)
                )
            )
        );
        volVectorField& mnf = tmnf.ref();
        mnf = constant::electromagnetic::mu0 * mag(H) * fvc::grad(mag(H));
        return tmnf;
    }
    */


    //NW 2017-11-13 21:11
    Info<< "Solving for U...\n" << endl;

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
//            ==
//            fvc::reconstruct
//            (
//                (
//                    mixture.surfaceTensionForce()
//                  - ghf*fvc::snGrad(rho)
//                  - fvc::snGrad(p_rgh + p_m)
//                  + constant::electromagnetic::mu0
//                  * fvc::interpolate(mag(H))*fvc::snGrad(mag(H))
//                ) * mesh.magSf()
//            )
        );

        //NW 2017-11-16 12:31
        Info<< "Correcting U...\n" << endl;
        fvOptions.correct(U);
    }
