scalar YagiAlphaBeta = 0.0;
scalar YagiPhi = 0.0;

if (kEffClosureModel_ == "Yagi")
{
    if
    (
        !porousMediaDict.found("YagiAlphaBeta")
     || !porousMediaDict.found("YagiPhi")
    )
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "ThermalConductivityClosureModel = Yagi requires "
            << "'YagiAlphaBeta' and 'YagiPhi'."
            << exit(FatalIOError);
    }

    YagiAlphaBeta = readScalar(porousMediaDict.lookup("YagiAlphaBeta"));
    YagiPhi = readScalar(porousMediaDict.lookup("YagiPhi"));
}

// Wakao mixing-conductivity coefficient
scalar Cmix = 0.1;

// Kunii stagnant-bed model parameters
scalar KuniiPhi = 0.0;
scalar KuniiGamma = 0.0;

if (kEffClosureModel_ == "Wakao")
{
    if (porousMediaDict.found("WakaoMixCoeff"))
    {
        Cmix = readScalar(porousMediaDict.lookup("WakaoMixCoeff"));
    }

    if (stagnantKeffClosureModel_ == "Kunii")
    {
        if
        (
            !porousMediaDict.found("KuniiPhi")
         || !porousMediaDict.found("KuniiGamma")
        )
        {
            FatalIOErrorInFunction(porousMediaDict)
                << "StagnantThermalConductivityClosureModel = Kunii "
                << "requires 'KuniiPhi' and 'KuniiGamma'."
                << exit(FatalIOError);
        }
        KuniiPhi   = readScalar(porousMediaDict.lookup("KuniiPhi"));
        KuniiGamma = readScalar(porousMediaDict.lookup("KuniiGamma"));
    }
}

// Create fields for prescribed porous-medium and solid-phase properties.
