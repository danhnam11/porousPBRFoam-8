const bool kEffClosureUsesSolidConductivity =
(
    kEffClosureModel_ == "Standard"
 || kEffClosureModel_ == "Yagi"
 || kEffClosureModel_ == "Wakao"
);

if (kEffClosureUsesSolidConductivity)
{
    if (!porousMediaDict.found("SolidThermalConductivity"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "The selected effective thermal-conductivity closure '"
            << kEffClosureModel_ << "' requires the entry "
            << "'SolidThermalConductivity' [W/m/K]."
            << exit(FatalIOError);
    }

    kappas =
        readScalar(porousMediaDict.lookup("SolidThermalConductivity"));
    alphas = kappas/Cps;
}

const bool kEffClosureUsesRe =
(
    kEffClosureModel_ == "Bunnell"
 || kEffClosureModel_ == "Yagi"
 || kEffClosureModel_ == "Wakao"
);

if (kEffClosureUsesRe)
{
    if (!porousMediaDict.found("SuperficialVelocity"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "The selected effective thermal-conductivity closure '"
            << kEffClosureModel_ << "' requires the entry "
            << "'SuperficialVelocity' [m/s] to evaluate the particle "
            << "Reynolds number."
            << exit(FatalIOError);
    }

    usup = readScalar(porousMediaDict.lookup("SuperficialVelocity"));
}
// Read and validate the stagnant-bed conductivity closure required by Wakao.
word stagnantKeffClosureModel_("None");

if (kEffClosureModel_ == "Wakao")
{
	if
    (
        !porousMediaDict.found
        (
            "StagnantThermalConductivityClosureModel"
        )
    )
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "ThermalConductivityClosureModel = Wakao requires entry "
            << "'StagnantThermalConductivityClosureModel'." << nl
            << "Valid options are: Krupiczka, Kunii, ZS"
            << exit(FatalIOError);
    }

    stagnantKeffClosureModel_ =
        word
        (
            porousMediaDict.lookup
            (
                "StagnantThermalConductivityClosureModel"
            )
        );
    if
    (
        stagnantKeffClosureModel_ != "Krupiczka"
     && stagnantKeffClosureModel_ != "Kunii"
     && stagnantKeffClosureModel_ != "ZS"
    )
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "Unknown StagnantThermalConductivityClosureModel: "
            << stagnantKeffClosureModel_ << nl
            << "Valid options are: Krupiczka, Kunii, ZS"
            << exit(FatalIOError);
    }
}

// Initialize the viscous and inertial resistance coefficients. The selected
// pressure-drop closure assigns their values below.
