#include "porousMedia/readProperties/readPorousMediaDict.H"

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
#include "porousMedia/readProperties/readPorousMediaProperties.H"

#include "porousMedia/readProperties/readThermalConductivityModel.H"

#include "porousMedia/readProperties/readThermalClosureInputs.H"

#include "porousMedia/createFields/createPorousMomentum.H"

#include "porousMedia/readProperties/readConstantConductivity.H"

#include "porousMedia/readProperties/readAnisotropicConductivity.H"

#include "porousMedia/readProperties/readThermalClosureCoeffs.H"

#include "porousMedia/createFields/createPorousMediaFields.H"

#include "porousMedia/createFields/createThermalConductivityTensor.H"

#include "porousMedia/createFields/createPorousTransportFields.H"

#include "porousMedia/createFields/createSpeciesDiffusivity.H"

autoPtr<psiReactionCatThermo> pCatThermo;
autoPtr<BasicSurfaceChemistryModel<psiReactionThermo, psiReactionCatThermo>> surfaceReaction;

#include "porousMedia/createFields/createPorousEnergyFields.H"

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
