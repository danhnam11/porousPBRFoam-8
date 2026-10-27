dimensionedScalar kEffRadial
(
    "kEffRadial",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    0.0
);

dimensionedScalar kEffAxial
(
    "kEffAxial",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    0.0
);

vector kEffAxialDirection(vector::zero);

if (kEffClosureModel_ == "ConstAnisotropic")
{
    if
    (
        !porousMediaDict.found("EffectiveRadialThermalConductivity")
     || !porousMediaDict.found("EffectiveAxialThermalConductivity")
     || !porousMediaDict.found("AxialDirection")
    )
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "ThermalConductivityClosureModel = ConstAnisotropic requires "
            << "'EffectiveRadialThermalConductivity', "
            << "'EffectiveAxialThermalConductivity' and 'AxialDirection'."
            << exit(FatalIOError);
    }

    kEffRadial = dimensionedScalar
    (
        "kEffRadial",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar
        (
            porousMediaDict.lookup("EffectiveRadialThermalConductivity")
        )
    );

    kEffAxial = dimensionedScalar
    (
        "kEffAxial",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar
        (
            porousMediaDict.lookup("EffectiveAxialThermalConductivity")
        )
    );

    kEffAxialDirection = vector(porousMediaDict.lookup("AxialDirection"));

    if (mag(kEffAxialDirection) < SMALL)
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "AxialDirection must be a non-zero vector."
            << exit(FatalIOError);
    }

    kEffAxialDirection /= mag(kEffAxialDirection);
}

// Yagi model parameters.
