// Update the porous-medium momentum resistance coefficient.
if (deltaPClosureModel_ == "None")
{
    porousMomentumCoeff =
        dimensionedScalar("zero", dimDensity/dimTime, 0.0);
}
else if (deltaPClosureModel_ == "Darcy")
{
    tmp<volScalarField> tmu = thermo.mu();
    const volScalarField& Intr_muF = tmu();

    porousMomentumCoeff =
        Intr_muF*DarcyCoeff;
}
else if
(
    deltaPClosureModel_ == "Forchheimer"
 || deltaPClosureModel_ == "Ergun"
)
{
    tmp<volScalarField> tmu = thermo.mu();
    const volScalarField& Intr_muF = tmu();

    porousMomentumCoeff =
        Intr_muF*DarcyCoeff
      + 0.5*Intr_rhoF*mag(Extr_UF)*ForchCoeff;
}

porousMomentumCoeff.correctBoundaryConditions();
