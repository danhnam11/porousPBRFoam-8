Info << "\n: create properties for porous media calculations " << endl;

// -----------------------------------------------------------------------------
// Read porous-media property dictionary
// -----------------------------------------------------------------------------

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

// default = off → energy equation solved
const Switch isothermal =
    porousMediaDict.lookupOrDefault<Switch>
    (
        "isothermal",
        false
    );

Info<< "Isothermal mode : "
    << (isothermal ? "ON" : "OFF")
    << nl << endl;

// porosity [-]
scalar por     (readScalar(porousMediaDict.lookup("Porosity")));
// particle diameter [m]
scalar dp      (readScalar(porousMediaDict.lookup("ParticleDiameter")));
// particle specific heat capacity [kg-m2/s2-K]
scalar Cps     (readScalar(porousMediaDict.lookup("SolidHeatCapacity")));
// Particle thermal conductivity [W/m-K] and diffusivity [m2/s]. They are
// loaded below only for closures that explicitly use the solid conductivity.
scalar kappas = 0.0;
scalar alphas = 0.0;
// particle density [kg/m3]
scalar rhos    (readScalar(porousMediaDict.lookup("SolidDensity")));
// Superficial velocity [m/s]. It is read below only by closure models that
// use the particle Reynolds number. Keep a zero initialization for models
// whose equations do not use Re.
scalar usup = 0.0;
// dummy parameter for PM mode
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

// Only these closures explicitly contain the solid thermal conductivity.
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

// Bunnell, Yagi and Wakao use Re = rho*usup*dp/mu. For the remaining
// effective-conductivity closures, SuperficialVelocity is not an input to the
// selected equation and therefore must not be required from the dictionary.
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

// Coefficient for mixing conductivity
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
    //dimensionedScalar("zero",dimless, 1.0),
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
    //dimensionedScalar("zero",dimless, 1.0),
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

volScalarField Pe
(
    IOobject
    (
        "Pe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

Info<< "Creating field effective mass diffusivity DiEff_i\n" << endl;
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

// ---------- Declarations only (no dictionary reads) ----------


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
