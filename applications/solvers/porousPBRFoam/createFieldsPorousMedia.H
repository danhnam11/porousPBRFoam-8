Info << "\n: create properties for porous media calculations " << endl;

// Read porous-medium model settings and material properties.

IOdictionary porousMediaDict
(
    IOobject
    (
        "porousMediaProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read and validate the solver mode.
if (!porousMediaDict.found("solverMode"))
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Missing required entry 'solverMode'." << nl
        << "Valid options are: PMSCHEM, PM"
        << exit(FatalIOError);
}

word solverMode_ = word(porousMediaDict.lookup("solverMode"));
Info<< "Selected solverMode = " << solverMode_ << endl;

// PMSCHEM: porous-media model with surface chemistry
// PM:      non-reacting porous-media model
if
( solverMode_ != "PMSCHEM" && solverMode_ != "PM" )
{   
    FatalErrorInFunction
        << "Invalid solverMode : " << solverMode_ << nl
        << "Available solverMode options are:" << nl
        << "    - PMSCHEM : porous media with surface chemistry" << nl 
        << "    - PM      : porous media without surface chemistry" << nl
        << exit(FatalError);
}

// Solve the energy equation unless isothermal operation is requested.
const Switch isothermal =
    porousMediaDict.lookupOrDefault<Switch>
    (
        "isothermal",
        false
    );

if (isothermal)
{
    Info<< "Isothermal mode: ON" << endl;
}

// Prescribed porous-medium and solid-phase properties.
// Porosity [-]
scalar por     (readScalar(porousMediaDict.lookup("Porosity")));
// Particle diameter [m]
scalar dp      (readScalar(porousMediaDict.lookup("ParticleDiameter")));
// Solid specific heat capacity [J kg^-1 K^-1]
scalar Cps     (readScalar(porousMediaDict.lookup("SolidHeatCapacity")));
// Solid thermal conductivity [W m^-1 K^-1]
scalar kappas  (readScalar(porousMediaDict.lookup("SolidThermalConductivity")));
// Solid enthalpy diffusivity, kappa_s/Cp_s [kg m^-1 s^-1]
scalar alphas = kappas/Cps;
// Solid density [kg m^-3]
scalar rhos    (readScalar(porousMediaDict.lookup("SolidDensity")));
// Reference superficial velocity used by closure correlations [m s^-1]
scalar usup    (readScalar(porousMediaDict.lookup("SuperficialVelocity")));
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
// Read and validate the effective-thermal-conductivity closure.
if (!porousMediaDict.found("ThermalConductivityClosureModel"))
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Missing required entry 'ThermalConductivityClosureModel'." << nl
        << "Valid options are: Const, Standard, Bunnell, Yagi, Wakao"
        << exit(FatalIOError);
}

word kEffClosureModel_
(
    porousMediaDict.lookup("ThermalConductivityClosureModel")
);

if
(
    kEffClosureModel_ != "Const"
 && kEffClosureModel_ != "Standard"
 && kEffClosureModel_ != "Bunnell"
 && kEffClosureModel_ != "Yagi"
 && kEffClosureModel_ != "Wakao"
)
{
    FatalIOErrorInFunction(porousMediaDict)
        << "Unknown ThermalConductivityClosureModel: "
        << kEffClosureModel_ << nl
        << "Valid options are: Const, Standard, Bunnell, Yagi, Wakao"
        << exit(FatalIOError);
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

// Coefficient multiplying the superficial velocity in the porous momentum
// source. Its local value is updated in updateForPorousMedia.H.
volScalarField porousMomentumCoeff
(
    IOobject
    (
        "porousMomentumCoeff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimDensity/dimTime, 0.0)
);

// Parameters used by the effective-thermal-conductivity closures.
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

volScalarField POR
(
    IOobject
    (
        "POR",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, por),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField FS2V
(
    IOobject
    (
        "FS2V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, fs2v),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField DP
(
    IOobject
    (
        "DP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 1, 0, 0, 0, 0, 0), dp),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField USUP
(
    IOobject
    (
        "USUP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 1, -1, 0, 0, 0, 0), usup),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField rhoS
(
    IOobject
    (
        "rhoS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -3, 0, 0, 0, 0, 0), rhos),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField CpS
(
    IOobject
    (
        "CpS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(0, 2, -2, -1, 0, 0, 0), Cps),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaS
(
    IOobject
    (
        "kappaS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), kappas),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alphaEffS
(
    IOobject
    (
        "alphaEffS",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), alphas),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField pAlphaEff
(
    IOobject
    (
        "pAlphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaEffStag
(
    IOobject
    (
        "kappaEffStag",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappaEffMix
(
    IOobject
    (
        "kappaEffMix",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, 1, -3, -1, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Re
(
    IOobject
    (
        "Re",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField Pr
(
    IOobject
    (
        "Pr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

PtrList<volScalarField> DiEff(Y.size());
forAll(DiEff, i)
{
    DiEff.set
    (
        i,
         new volScalarField
         (
             IOobject
             (
               "DiEff_"+composition.species()[i],
               runTime.timeName(),
               mesh,
               IOobject::NO_READ,
               IOobject::AUTO_WRITE
             ),
	     mesh,
	     dimensionedScalar("DiEff"+composition.species()[i], dimensionSet(0, 2, -1, 0, 0), 0.)
         )
    );
}

// Declare surface-chemistry objects and energy-source fields. Their model-
// specific initialization is performed only when PMSCHEM mode is selected.

autoPtr<psiReactionCatThermo> pCatThermo;
autoPtr<BasicSurfaceChemistryModel<psiReactionThermo, psiReactionCatThermo>> surfaceReaction;

volScalarField heGas
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

Info<< nl
    << "Porous-media closure models:" << nl
    << "    Pressure drop                     : "
    << deltaPClosureModel_ << nl
    << "    Effective thermal conductivity    : "
    << kEffClosureModel_ << nl;

if (kEffClosureModel_ == "Wakao")
{
    Info<< "    Stagnant-bed thermal conductivity : "
        << stagnantKeffClosureModel_ << nl;
}

Info<< "    Effective species diffusivity     : "
    << DiEffClosureModel_ << endl;
