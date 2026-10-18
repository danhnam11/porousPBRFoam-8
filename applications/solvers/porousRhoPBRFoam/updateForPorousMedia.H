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

if (!isothermal)
{
// Update the particle Reynolds and fluid Prandtl numbers used by thermal
// closure models.
Re = (DP*USUP)*rho/mu;
Pr = Cp*mu/kappa;

// Initialize the effective thermal-transport fields using fluid properties.
kappaEff = kappa;
pAlphaEff = thermophysicalTransport->alphaEff();

if (kEffClosureModel_ == "Const")
{
    kappaEff  = kEffConst;
    pAlphaEff = kEffConst/Cp;
}
else if (kEffClosureModel_ == "ConstAnisotropic")
{
    // kappaEff remains the wall-normal (radial) conductivity used by the
    // porousMediaHeatTransferWall boundary condition.
    kappaEff  = kEffRadial;
    pAlphaEff = kEffRadial/Cp;
}
else if (kEffClosureModel_ == "Standard")
{
    alphaEff  = thermophysicalTransport->alphaEff();
    pAlphaEff = por*alphaEff + (1-por)*alphaEffS*CpS/Cp;
    kappaEff  = por*kappa + (1-por)*kappaS;
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
    // ratios
    const volScalarField Rsk = kappaS/kappa; // ks/kf
    const volScalarField Rfk = kappa/kappaS; // kf/ks

    kappaEffStag = kappa; // fallback

    if (stagnantKeffClosureModel_ == "Krupiczka")
    {
        // a = 0.280 - 0.757*log(por) - 0.057*log(ks/kf)
	volScalarField logPOR = Foam::log10(POR);
        volScalarField logR   = Foam::log10(Rsk);


	volScalarField a = 0.280 - 0.757*logPOR - 0.057*logR;

        // keff0/kf = (ks/kf)^a  => keff0 = kf * (ks/kf)^a
        kappaEffStag = kappa * Foam::pow(Rsk, a);
    }
    else if (stagnantKeffClosureModel_ == "KuniiSmith")
    {
        volScalarField denom = KuniiPhi + KuniiGamma*Rfk;
        // keff0/kf = por + (1-eps)/(phi + gamma*(kf/ks))
        denom = Foam::max(denom, dimensionedScalar("small", denom.dimensions(), SMALL));
        //
        volScalarField keff0_over_kf = POR + (1.0 - POR)/denom;

        kappaEffStag = kappa * keff0_over_kf;

    }
    else if (stagnantKeffClosureModel_ == "ZS")
    {
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
    else
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "Unknown stagnantKeffModel: " << stagnantKeffClosureModel_ << nl
            << "Valid: Krupiczka, KuniiSmith, ZS"
            << exit(FatalIOError);
    }

        kappaEffMix = Cmix * kappa * Re * Pr;

        kappaEff  = kappaEffStag + kappaEffMix;
        pAlphaEff = kappaEff/Cp;
}

else
{
    FatalErrorInFunction
        << "Unknown kEffMode : " << kEffClosureModel_ << nl
        << "Available options are:" << nl
        << "  1. Const" << nl
        << "  2. ConstAnisotropic" << nl
        << "  3. Standard" << nl
        << "  4. Bunnell" << nl
        << "  5. Yagi" << nl
        << "  6. Wakao" << nl
        << exit(FatalError);
}

// Convert the scalar closures to isotropic tensors.  For Wen's 2DADPF
// closure, replace the axial component with the independently prescribed
// effective axial conductivity.
pAlphaEffTensor = pAlphaEff*symmTensor::I;

if (kEffClosureModel_ == "ConstAnisotropic")
{
    pAlphaEffTensor =
        (kEffRadial/Cp)*symmTensor::I
      + ((kEffAxial - kEffRadial)/Cp)
       *sqr(kEffAxialDirection);
}
}

forAll(DiEff, i)
{
    DiEff[i] =
        sqrt((POR+SMALL))*(thermophysicalTransport->alphaEff()/rho)
      + (mag(U)*DP)/(12.0*(POR+SMALL));

}
