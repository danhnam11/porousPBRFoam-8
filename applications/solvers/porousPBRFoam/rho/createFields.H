#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
rhoReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Intr_YF = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

volScalarField Intr_rhoF
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField Extr_UF
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField Intr_muF
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);

volScalarField Intr_CpF
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);

volScalarField Intr_kappaF
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);

volScalarField alphaEff
(
    IOobject
    (
        "alphaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimensionSet(1, -1, -1, 0, 0, 0, 0), 0.0)
);


volScalarField& Intr_pF = thermo.p();

// Compatibility references for standard OpenFOAM headers.
// These names refer to the same objects; no copying or averaging.
volVectorField& U = Extr_UF;
volScalarField& rho = Intr_rhoF;
#include "compressibleCreatePhi.H"

pressureControl pressureControl(Intr_pF, Intr_rhoF, pimple.dict(), false);

mesh.setFluxRequired(Intr_pF.name());


Info << "Creating turbulence model.\n" << nl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        Intr_rhoF,
        Extr_UF,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<rhoReactionThermophysicalTransportModel>
thermophysicalTransport
(
    rhoReactionThermophysicalTransportModel::New
    (
        turbulence(),
        thermo
    )
);

Info<< "Creating reaction model\n" << endl;
autoPtr<CombustionModel<rhoReactionThermo>> reaction
(
    CombustionModel<rhoReactionThermo>::New(thermo, turbulence())
);


Info<< "Creating field dpdt\n" << endl;
volScalarField Intr_dpdtF
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar(Intr_pF.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(Extr_UF));

dimensionedScalar initialMass = fvc::domainIntegrate(Intr_rhoF);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Intr_YF, i)
{
    fields.add(Intr_YF[i]);
}
fields.add(thermo.he());

#include "createMRF.H"
#include "createFvOptions.H"
