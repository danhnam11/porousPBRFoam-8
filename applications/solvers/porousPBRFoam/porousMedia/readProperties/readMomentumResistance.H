// Momentum-exchange source: select the closure and read its coefficients.
validatePorousMomentumInputs(porousMediaDict);
word deltaPClosureModel_(porousMediaDict.lookup("pressureDropClosureModel"));

dimensionedScalar DarcyCoeff
(
    "DarcyCoeff",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),  // [1/m2]
    0.0
);

dimensionedScalar ForchCoeff
(
    "ForchCoeff",
    dimensionSet(0, -1, 0, 0, 0, 0, 0),  // [1/m]
    0.0
);

// Read or evaluate the resistance coefficients for the selected closure.
if (deltaPClosureModel_ == "None")
{
    Info<< "No porous pressure-drop source is applied." << endl;
}
else if (deltaPClosureModel_ == "Darcy")
{
    if (!porousMediaDict.found("DarcyCoeff"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "pressureDropClosureModel = Darcy requires entry "
            << "'DarcyCoeff'."
            << exit(FatalIOError);
    }

    DarcyCoeff = dimensionedScalar
    (
        "DarcyCoeff",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        readScalar(porousMediaDict.lookup("DarcyCoeff"))
    );

    ForchCoeff = dimensionedScalar
    (
        "ForchCoeff",
        dimensionSet(0, -1, 0, 0, 0, 0, 0),
        0.0
    );

    Info<< "DarcyCoeff = " << DarcyCoeff.value() << " [1/m2]" << endl;
}
else if (deltaPClosureModel_ == "Forchheimer")
{
    if (!porousMediaDict.found("DarcyCoeff"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "pressureDropClosureModel = Forchheimer requires entry "
            << "'DarcyCoeff'."
            << exit(FatalIOError);
    }

    if (!porousMediaDict.found("ForchCoeff"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "pressureDropClosureModel = Forchheimer requires entry "
            << "'ForchCoeff'."
            << exit(FatalIOError);
    }

    DarcyCoeff = dimensionedScalar
    (
        "DarcyCoeff",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        readScalar(porousMediaDict.lookup("DarcyCoeff"))
    );

    ForchCoeff = dimensionedScalar
    (
        "ForchCoeff",
        dimensionSet(0, -1, 0, 0, 0, 0, 0),
        readScalar(porousMediaDict.lookup("ForchCoeff"))
    );

    Info<< "DarcyCoeff = " << DarcyCoeff.value() << " [1/m2]" << endl;
    Info<< "ForchCoeff = " << ForchCoeff.value() << " [1/m]" << endl;
}
else if (deltaPClosureModel_ == "Ergun")
{
    // The Ergun coefficients depend only on porosity and particle diameter and
    // therefore need to be evaluated only during field initialization.
    const dimensionedScalar particleDiameter
    (
        "ParticleDiameter",
        dimLength,
        dp
    );

    DarcyCoeff = dimensionedScalar
    (
        "DarcyCoeff",
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        150.0*sqr(1.0 - por)
       /(pow3(por)*sqr(particleDiameter.value()))
    );

    ForchCoeff = dimensionedScalar
    (
        "ForchCoeff",
        dimensionSet(0, -1, 0, 0, 0, 0, 0),
        3.5*(1.0 - por)
       /(pow3(por)*particleDiameter.value())
    );

    Info<< "Ergun model:" << nl
        << "  DarcyCoeff       = " << DarcyCoeff.value() << " [1/m2]" << nl
        << "  ForchCoeff       = " << ForchCoeff.value() << " [1/m]" << endl;
}

