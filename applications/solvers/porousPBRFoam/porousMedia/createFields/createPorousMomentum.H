// Coefficient multiplying the superficial velocity in the porous momentum
// source. Its local value is updated in updateForPorousMedia.H.
volScalarField porousMomentumCoeff
(
    IOobject
    (
        "porousMomentumCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);
