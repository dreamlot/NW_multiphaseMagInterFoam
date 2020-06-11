
#include "electromagneticConstants.H"


Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//NW20170920
Info<< "Reading field psiMag\n" << endl;
volScalarField psiMag
(
    IOobject
    (
        "psiMag",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



#include "createPhi.H"

NW_multiphaseMagMixture mixture(U, phi, psiMag);

// Need to store rho for ddt(rho, U)
Info<< "Reading field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.rho()
);
rho.oldTime();

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture)
);



Info<< "Reading field muMag\n" << endl;
volScalarField muMag
(
    IOobject
    (
        "muMag",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.muMag()
);
muMag.oldTime();

//NW 2017-11-14 21:38
Info<< "Reading solid phase flag\n" << endl;
volScalarField isSolid
(
    IOobject
    (
        "isSolid",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mixture.isSolid()
);
isSolid.oldTime();


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

// output: test codes NW20171017



volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"

/*
volVectorField mnfMag
(
    IOobject
    (
        "mnfMag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("mnfMag", dimensionSet(1,-5,-2,0,0,0,0), vector(0.0,0.0,0.0))


);
*/



Info<< "Reading field H\n" << endl;
volVectorField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //mesh
    fvc::reconstruct(fvc::snGrad(psiMag)*mesh.magSf())
);

Info<< "Reading field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //mesh
    constant::electromagnetic::mu0
    *muMag*fvc::reconstruct(fvc::snGrad(psiMag)*mesh.magSf())
);

// fluid-magnetic pressure
volScalarField p_m
(
    IOobject
    (
        "p_m",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //p_rgh
    //0.5*constant::electromagnetic::mu0 * (muMag-1.0) * H & H
    0.5*constant::electromagnetic::mu0 * (muMag-1.0) * magSqr(H)
);

volScalarField magHMag
(
    IOobject
    (
        "magHMag",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mag(H)
);
/*
volScalarField fnMag
(
    IOobject
    (
        "fnMag",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    //constant::electromagnetic::mu0 * mag(H)*fvc::snGrad(mag(H)) * mesh.magSf()
    constant::electromagnetic::mu0
    *fvc::snGrad(fvc::reconstruct(mag(fvc::snGrad(psiMag)*mesh.magSf())))
);
*/
