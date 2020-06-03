
while (pimple.correctNonOrthogonal())
{
    solve(fvm::laplacian(muMag, psiMag));			//ningyu 20161027
}
H = fvc::reconstruct(fvc::snGrad(psiMag)*mesh.magSf());
B = constant::electromagnetic::mu0*muMag*H;

// magnetic normal pressure
volVectorField MMag
(
    "MMag",
    (muMag-1.0)*H
);
volScalarField magMMag = mag(MMag);
magHMag = mag(H);
//surfaceScalarField gradHmag = fvc::grad(magHMag);
//surfaceScalarField mnfMag = magHMag*gradHmag;
//volVectorField mnfMag = MMag & fvc::grad(H);

/*
tmp<surfaceVectorField> tmnfMag
(
    new surfaceVectorField
    (
        IOobject
        (
            "magneticNormalForce",
            mesh.time().timeName(),
            mesh
        )
        ,mesh,
        dimensionedVector("magneticNormalForce",dimensionSet(1,-5,-2,0,0,0,0), vector(0.0,0.0,0.0))
    )
);
surfaceVectorField& mnf = tmnfMag.ref();
//mnf = dimensionedScalar("tmpName",dimless,1.0) * linearInterpolate(mnfMag);
mnf = fvc::interpolate(mnfMag);
*/
/*
tmp<surfaceScalarField> tmnfMag
(
    new surfaceScalarField
    (
        IOobject
        (
            "magneticNormalForce",
            mesh.time().timeName(),
            mesh
        )
        ,mesh,
        dimensionedScalar("magneticNormalForce",dimensionSet(1,-5,-2,0,0,0,0), 0.0)
    )
);
surfaceScalarField& mnf = tmnfMag.ref();
//mnf = dimensionedScalar("tmpName",dimless,1.0) * linearInterpolate(mnfMag);
//mnf = fvc::interpolate(mag(mnfMag));
mnf =  magHMag*fvc::snGrad(magHMag);
*/

//volScalarField fn =  magHMag*fvc::snGrad(magHMag);

// Fluid-magnetic pressure. The following two expressions are equivalent.
//p_m = 0.5*constant::electromagnetic::mu0 * (muMag-1.0) * magHMag * magHMag;
p_m = 0.5*constant::electromagnetic::mu0 * (muMag-1.0) * H & H;
