Info<< "Reading catalyst thermophysical properties\n" << endl;

// Create the catalyst thermophysical model.
pCatThermo.reset
(
    psiReactionCatThermo::New(mesh).ptr()
);

psiReactionCatThermo& catThermo = pCatThermo();
catThermo.validate(args.executable(), "h", "e");

Info<< "Catalyst thermophysical properties created\n" << endl;
Info<< "Creating catalyst surface-chemistry model\n" << endl;

// Create the heterogeneous surface-chemistry model.
surfaceReaction =
    BasicSurfaceChemistryModel<psiReactionThermo, psiReactionCatThermo>::New
    (
        thermo,
        catThermo
    );

// Read the surface-chemistry properties.
IOdictionary surfaceChemistryDict
(
    IOobject
    (
        "surfaceChemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Read the catalytic surface-to-volume ratio [m^-1].
fs2v =
    readScalar
    (
        surfaceChemistryDict.lookup("SurfaceToVolumeRatio")
    );

Info << "SurfaceToVolumeRatio = " << fs2v << " [1/m]" << endl;

// Assign the surface-to-volume ratio to the computational field.
FS2V.primitiveFieldRef() = fs2v;
FS2V.correctBoundaryConditions();

Info<< "Catalyst surface-chemistry model created\n" << endl;
