volScalarField POR
(
    IOobject
    (
        "POR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, por),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField FS2V
(
    IOobject
    (
        "FS2V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, fs2v),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField DP
(
    IOobject
    (
        "DP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 1, 0, 0, 0, 0, 0), dp),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField USUP
(
    IOobject
    (
        "USUP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 1, -1, 0, 0, 0, 0), usup),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Intr_rhoS
(
    IOobject
    (
        "rhoS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -3, 0, 0, 0, 0, 0), rhos),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Intr_CpS
(
    IOobject
    (
        "CpS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 2, -2, -1, 0, 0, 0), Cps),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Intr_kappaS
(
    IOobject
    (
        "kappaS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), kappas),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alphaEffS
(
    IOobject
    (
        "alphaEffS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), alphas),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField pAlphaEff
(
    IOobject
    (
        "pAlphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);
