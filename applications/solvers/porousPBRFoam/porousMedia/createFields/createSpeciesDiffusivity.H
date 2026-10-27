PtrList<volScalarField> DiEff(Intr_YF.size());
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
