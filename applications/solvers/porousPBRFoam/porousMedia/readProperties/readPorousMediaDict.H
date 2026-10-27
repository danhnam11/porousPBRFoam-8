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
