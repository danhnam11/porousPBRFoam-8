volSymmTensorField pAlphaEffTensor
(
    IOobject
    (
        "pAlphaEffTensor",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedSymmTensor
    (
        "zero",
        dimensionSet(1, -1, -1, 0, 0, 0, 0),
        symmTensor::zero
    )
);
