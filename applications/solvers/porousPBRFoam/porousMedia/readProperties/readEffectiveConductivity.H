// Effective thermal conductivity: selection, input validation and coefficients.
// Keep declarations in the enclosing solver scope for subsequent field creation.
word kEffClosureModel_("None");
word stagnantKeffClosureModel_("None");
scalar kappas = 0.0;
scalar alphas = 0.0;
scalar usup = 0.0;
scalar YagiAlphaBeta = 0.0;
scalar YagiPhi = 0.0;
scalar Cmix = 0.1;
scalar KuniiPhi = 0.0;
scalar KuniiGamma = 0.0;

// Thermal conductivity dimensions: [W/m/K].
dimensionedScalar kEffConst
(
    "kEffConst", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0
);
dimensionedScalar kEffRadial
(
    "kEffRadial", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0
);
dimensionedScalar kEffAxial
(
    "kEffAxial", dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0
);
vector kEffAxialDirection(vector::zero);

if (!isothermal)
{
    // Print the selected model's required keywords before any coefficient lookup.
    validatePorousThermalInputs(porousMediaDict);
    kEffClosureModel_ = word(porousMediaDict.lookup("ThermalConductivityClosureModel"));

    if
    (
        kEffClosureModel_ == "Standard"
     || kEffClosureModel_ == "Yagi"
     || kEffClosureModel_ == "Wakao"
    )
    {
        kappas = readScalar(porousMediaDict.lookup("SolidThermalConductivity"));
        alphas = kappas/Cps;
    }
    if
    (
        kEffClosureModel_ == "Bunnell"
     || kEffClosureModel_ == "Yagi"
     || kEffClosureModel_ == "Wakao"
    )
    {
        usup = readScalar(porousMediaDict.lookup("SuperficialVelocity"));
    }

    if (kEffClosureModel_ == "Const")
    {
        kEffConst.value() = readScalar(porousMediaDict.lookup("EffectiveThermalConductivity"));
    }
    else if (kEffClosureModel_ == "ConstAnisotropic")
    {
        kEffRadial.value() = readScalar(porousMediaDict.lookup("EffectiveRadialThermalConductivity"));
        kEffAxial.value() = readScalar(porousMediaDict.lookup("EffectiveAxialThermalConductivity"));
        kEffAxialDirection = vector(porousMediaDict.lookup("AxialDirection"));
        if
        (
            !std::isfinite(kEffAxialDirection.x())
         || !std::isfinite(kEffAxialDirection.y())
         || !std::isfinite(kEffAxialDirection.z())
         || !std::isfinite(mag(kEffAxialDirection))
         || mag(kEffAxialDirection) < SMALL
        )
        {
            FatalIOErrorInFunction(porousMediaDict)
                << "ThermalConductivityClosureModel = ConstAnisotropic requires "
                << "AxialDirection (e_z) to be a finite non-zero vector."
                << exit(FatalIOError);
        }
        kEffAxialDirection /= mag(kEffAxialDirection);
    }
    else if (kEffClosureModel_ == "Yagi")
    {
        YagiAlphaBeta = readScalar(porousMediaDict.lookup("YagiAlphaBeta"));
        YagiPhi = readScalar(porousMediaDict.lookup("YagiPhi"));
    }
    else if (kEffClosureModel_ == "Wakao")
    {
        stagnantKeffClosureModel_ = word
        (
            porousMediaDict.lookup("StagnantThermalConductivityClosureModel")
        );
        Cmix = porousMediaDict.lookupOrDefault<scalar>("WakaoMixCoeff", 0.1);
        if (stagnantKeffClosureModel_ == "Kunii")
        {
            KuniiPhi = readScalar(porousMediaDict.lookup("KuniiPhi"));
            KuniiGamma = readScalar(porousMediaDict.lookup("KuniiGamma"));
        }
    }
}
else
{
    Info<< "Isothermal mode: thermal-conductivity closure input is skipped."
        << nl << endl;
}
