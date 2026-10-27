scalar por     (readScalar(porousMediaDict.lookup("Porosity")));
// Particle diameter [m]
scalar dp      (readScalar(porousMediaDict.lookup("ParticleDiameter")));
// Solid specific heat capacity [J kg^-1 K^-1]
scalar Cps     (readScalar(porousMediaDict.lookup("SolidHeatCapacity")));
// Solid thermal conductivity and enthalpy diffusivity. They are loaded below
// only for closures that explicitly use the solid conductivity.
scalar kappas = 0.0;
scalar alphas = 0.0;
// Solid density [kg m^-3]
scalar rhos    (readScalar(porousMediaDict.lookup("SolidDensity")));
// Reference superficial velocity used by Reynolds-based closure correlations.
// It is loaded below only when the selected closure requires it.
scalar usup = 0.0;
// Catalytic surface-to-volume ratio remains zero in non-reacting PM mode.
scalar fs2v = 0.0;

// Read and validate the pressure-drop closure.
if (!porousMediaDict.found("pressureDropClosureModel"))
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Missing required entry 'pressureDropClosureModel'." << nl
        << "Valid options are: None, Darcy, Forchheimer, Ergun"
        << exit(FatalIOError);
}

word deltaPClosureModel_
(
    porousMediaDict.lookup("pressureDropClosureModel")
);

if
(
    deltaPClosureModel_ != "None"
 && deltaPClosureModel_ != "Darcy"
 && deltaPClosureModel_ != "Forchheimer"
 && deltaPClosureModel_ != "Ergun"
)
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Unknown pressureDropClosureModel: "
        << deltaPClosureModel_ << nl
        << "Valid options are: None, Darcy, Forchheimer, Ergun"
        << exit(FatalIOError);
}

// Read and validate the effective-species-diffusivity closure.
if (!porousMediaDict.found("DiffusivityClosureModel"))
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Missing required entry 'DiffusivityClosureModel'." << nl
        << "Valid option is: BruggemanWilhelm"
        << exit(FatalIOError);
}

word DiEffClosureModel_
(
    porousMediaDict.lookup("DiffusivityClosureModel")
);

if (DiEffClosureModel_ != "BruggemanWilhelm")
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Unknown DiffusivityClosureModel: "
        << DiEffClosureModel_ << nl
        << "Valid option is: BruggemanWilhelm"
        << exit(FatalIOError);
}
// Thermal closure data are meaningful only when the energy equation is
// active. "None" is an internal sentinel and is not a user-selectable model.
word kEffClosureModel_("None");
