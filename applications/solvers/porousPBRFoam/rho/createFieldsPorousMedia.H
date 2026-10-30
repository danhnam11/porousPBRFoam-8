#include "porousMedia/readProperties/readPorousMediaDict.H"

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
#include "porousMedia/readProperties/readPorousMediaProperties.H"
#include "porousMedia/readProperties/readMomentumResistance.H"
#include "porousMedia/readProperties/readSpeciesDiffusivity.H"
#include "porousMedia/readProperties/readEffectiveConductivity.H"

#include "porousMedia/createFields/createPorousMomentum.H"

#include "porousMedia/createFields/createPorousMediaFields.H"

#include "porousMedia/createFields/createThermalConductivityTensor.H"

#include "porousMedia/createFields/createPorousTransportFields.H"

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
#include "porousMedia/createFields/createSpeciesDiffusivity.H"

#include "porousMedia/createFields/createPorousEnergyFields.H"

