dimensionedScalar kEffConst
(
    "kEffConst",
    dimensionSet(1, 1, -3, -1, 0, 0, 0),
    0.0
);

if (kEffClosureModel_ == "Const")
{
    if (!porousMediaDict.found("EffectiveThermalConductivity"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "kEffModel = Const requires entry "
            << "'EffectiveThermalConductivity'."
            << exit(FatalIOError);
    }

    kEffConst = dimensionedScalar
    (
        "kEffConst",
        dimensionSet(1, 1, -3, -1, 0, 0, 0),
        readScalar(porousMediaDict.lookup("EffectiveThermalConductivity"))
    );
}

// Yagi model parameters.
