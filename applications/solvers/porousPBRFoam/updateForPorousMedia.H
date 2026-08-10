// Update the porous-medium momentum resistance coefficient.
if (deltaPClosureModel_ == "None")
{
    porousMomentumCoeff =
        dimensionedScalar("zero", dimDensity/dimTime, 0.0);
}
else if (deltaPClosureModel_ == "Darcy")
{
    tmp<volScalarField> tmu = thermo.mu();
    const volScalarField& mu = tmu();

    porousMomentumCoeff =
        mu*DarcyCoeff;
}
else if
(
    deltaPClosureModel_ == "Forchheimer"
 || deltaPClosureModel_ == "Ergun"
)
{
    tmp<volScalarField> tmu = thermo.mu();
    const volScalarField& mu = tmu();

    porousMomentumCoeff =
        mu*DarcyCoeff
      + 0.5*rho*mag(U)*ForchCoeff;
}

porousMomentumCoeff.correctBoundaryConditions();

// Update the particle Reynolds and fluid Prandtl numbers.
Re = (DP*USUP)*rho/mu;
Pr = Cp*mu/kappa;

// Initialize the effective thermal-transport fields using fluid properties.
kappaEff = kappa;
pAlphaEff = thermophysicalTransport->alphaEff();

// Update the effective thermal conductivity using the selected closure.
if (kEffClosureModel_ == "Const")
{
    kappaEff  = kEffConst;
    pAlphaEff = kEffConst/Cp;
}
else if (kEffClosureModel_ == "Standard")
{
    kappaEff  = por*kappa + (1-por)*kappaS;
    pAlphaEff = kappaEff/Cp;
}
else if (kEffClosureModel_ == "Bunnell")
{
    kappaEff  = kappa*(5+0.061*Re);
    pAlphaEff = kappaEff/Cp;
}
else if (kEffClosureModel_ == "Yagi")
{
    kappaEff  = kappa*((1.0 - POR)*kappaS/(kappa + YagiPhi*kappaS) + (YagiAlphaBeta*Re*Pr));
    pAlphaEff = kappaEff/Cp;
}
else if (kEffClosureModel_ == "Wakao")
{
    // Ratios of solid-to-fluid and fluid-to-solid thermal conductivity.
    const volScalarField Rsk = kappaS/kappa; // k_s/k_f
    const volScalarField Rfk = kappa/kappaS; // k_f/k_s

    // Initialize the stagnant contribution with the fluid conductivity.
    kappaEffStag = kappa;

    if (stagnantKeffClosureModel_ == "Krupiczka")
    {
	// Evaluate the Krupiczka correlation exponent.
	volScalarField logPOR = Foam::log10(POR);
        volScalarField logR   = Foam::log10(Rsk);

	volScalarField a = 0.280 - 0.757*logPOR - 0.057*logR;

        // Evaluate k_eff^0/k_f = (k_s/k_f)^a.
        kappaEffStag = kappa * Foam::pow(Rsk, a);
    }
    else if (stagnantKeffClosureModel_ == "Kunii")
    {
        // Evaluate the denominator in the Kunii correlation.
        volScalarField denom = KuniiPhi + KuniiGamma*Rfk;

        // Prevent division by values smaller than the numerical tolerance.
        denom = Foam::max(denom, dimensionedScalar("small", denom.dimensions(), SMALL));

        // Evaluate k_eff^0/k_f = eps + (1-eps)/(phi + gamma*k_f/k_s).
        volScalarField keff0_over_kf = POR + (1.0 - POR)/denom;

        kappaEffStag = kappa * keff0_over_kf;
    }
    else if (stagnantKeffClosureModel_ == "ZS")
    {
        // Evaluate the intermediate terms of the Zehner--Schlunder model.
        const volScalarField A   = kappa/kappaS;
        const volScalarField B   = 1.25*Foam::pow((1.0 - POR)/POR, 10.0/9.0);

        const volScalarField s   = Foam::sqrt(1.0 - POR);
        volScalarField den = 1.0 - A*B;

        den = Foam::max
        (
            den,
            dimensionedScalar("small", den.dimensions(), SMALL)
        );

        volScalarField T1 = 1.0 - s;

        volScalarField T2 = 2.0*s/den;

        volScalarField T3 = (1.0 - A)*B/(den*den) * Foam::log(1.0/(A*B));

        volScalarField T4 = -(B + 1.0)/2.0;

        volScalarField T5 = -(B - 1.0)/den;

        volScalarField factor = T1 + T2*(T3 + T4 + T5);

        kappaEffStag = kappa * factor;
    }

    // Add the flow-dependent contribution of the Wakao correlation.
    kappaEffMix = Cmix * kappa * Re * Pr;

    kappaEff  = kappaEffStag + kappaEffMix;
    pAlphaEff = kappaEff/Cp;
}

// Update the effective species diffusivities using the selected closure.
if (DiEffClosureModel_ == "BruggemanWilhelm")
{
    forAll(DiEff, i)
    {
        DiEff[i] =
            sqrt(POR + SMALL)
           *(thermophysicalTransport->alphaEff()/rho)
          + (mag(U)*DP)/(12.0*(POR + SMALL));
    }
}