volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaEffStag
(
    IOobject
    (
        "kappaEffStag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaEffMix
(
    IOobject
    (
        "kappaEffMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Re
(
    IOobject
    (
        "Re",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Pr
(
    IOobject
    (
        "Pr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);
