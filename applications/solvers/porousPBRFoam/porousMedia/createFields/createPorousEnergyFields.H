volScalarField Intr_heGasF
(   
    IOobject
    (
        "heGas",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.he()
);

volScalarField QdotSurface
(
    IOobject
    (
        "QdotSurface",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
);
