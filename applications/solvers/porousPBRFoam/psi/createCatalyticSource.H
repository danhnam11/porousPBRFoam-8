// Create the catalytic species-source interface for psi thermodynamics.
// PM returns a zero matrix; PMSCHEM retains the existing source matrix.
auto catalyticSpeciesSource =
    [&](volScalarField& Intr_YiF) -> tmp<fvScalarMatrix>
    {
        tmp<fvScalarMatrix> tRcat
        (
            new fvScalarMatrix(Intr_YiF, dimMass/dimTime)
        );

        if (surfaceReaction.valid())
        {
            tRcat.ref() += FS2V*surfaceReaction->RFromSCHEM(Intr_YiF);
        }

        return tRcat;
    };
