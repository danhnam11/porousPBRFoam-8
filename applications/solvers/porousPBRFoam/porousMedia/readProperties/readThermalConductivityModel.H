if (!isothermal)
{
    if (!porousMediaDict.found("ThermalConductivityClosureModel"))
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "Missing required entry 'ThermalConductivityClosureModel'." << nl
            << "Valid options are: Const, ConstAnisotropic, Standard, Bunnell, Yagi, Wakao"
            << exit(FatalIOError);
    }

    kEffClosureModel_ =
        word(porousMediaDict.lookup("ThermalConductivityClosureModel"));

    if
    (
        kEffClosureModel_ != "Const"
     && kEffClosureModel_ != "ConstAnisotropic"
     && kEffClosureModel_ != "Standard"
     && kEffClosureModel_ != "Bunnell"
     && kEffClosureModel_ != "Yagi"
     && kEffClosureModel_ != "Wakao"
    )
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "Unknown ThermalConductivityClosureModel: "
            << kEffClosureModel_ << nl
            << "Valid options are: Const, ConstAnisotropic, Standard, Bunnell, Yagi, Wakao"
            << exit(FatalIOError);
    }
}
else
{
    Info<< "Isothermal mode: thermal-conductivity closure input is skipped."
        << nl << endl;
}
