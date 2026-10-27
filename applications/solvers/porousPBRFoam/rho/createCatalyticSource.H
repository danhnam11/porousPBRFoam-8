// The rho path supports PM and supplies no catalytic species source.
auto catalyticSpeciesSource =
    [](volScalarField& Intr_YiF) -> tmp<fvScalarMatrix>
    {
        return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix(Intr_YiF, dimMass/dimTime)
        );
    };
