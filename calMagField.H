if (!args.optionFound("noH") || args.optionFound("HdotGradH"))
    {
        volVectorField H
        (
            IOobject
            (
                "H",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            -fvc::reconstruct(fvc::snGrad(psiMag)*mesh.magSf())
        );
/*
        if (!args.optionFound("noH"))
        {
            Info<< nl
                << "Creating field H for time "
                << runTime.timeName() << endl;

            H.write();
        }
*/

        if (args.optionFound("HdotGradH"))
        {
            Info<< nl
                << "Creating field HdotGradH for time "
                << runTime.timeName() << endl;

            volVectorField HdotGradH
            (
                IOobject
                (
                    "HdotGradH",
                    runTime.timeName(),
                    mesh
                ),
                H & fvc::grad(H)
            );

//            HdotGradH.write();
        }
    }

    if (!args.optionFound("noB"))
    {
        Info<< nl
            << "Creating field B for time "
            << runTime.timeName() << endl;

        volVectorField B
        (
            IOobject
            (
                "B",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            constant::electromagnetic::mu0
            *muMag*fvc::reconstruct(fvc::snGrad(psiMag)*mesh.magSf())
        );

        //B.write();
    }


/*
    if (!args.optionFound("noB"))
    {
        Info<< nl
            << "Creating field snGrad_psi for time "
            << runTime.timeName() << endl;

        volVectorField tmpB
        (
            IOobject
            (
                "tmpB",
                runTime.timeName(),
                mesh
            ),
			fvc::reconstruct(fvc::snGrad(psi))
//            constant::electromagnetic::mu0
//           *fvc::reconstruct(murf*fvc::snGrad(psi)*mesh.magSf() + murf*Mrf)
        );

		tmpB.write();
	}

*/

/*
    if (!args.optionFound("noB"))
    {
        Info<< nl
            << "Creating field mesh.magSf for time "
            << runTime.timeName() << endl;

        volVectorField tmpB2
        (
            IOobject
            (
                "tmpB2",
                runTime.timeName(),
                mesh
            ),
			fvc::reconstruct(mesh.magSf())
//            constant::electromagnetic::mu0
//           *fvc::reconstruct(murf*fvc::snGrad(psi)*mesh.magSf() + murf*Mrf)
        );

        tmpB2.write();
    }
*/
