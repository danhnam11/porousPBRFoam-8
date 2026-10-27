if (!isothermal)
{
// Update the particle Reynolds and fluid Prandtl numbers used by thermal
// closure models.
#include "porousMedia/updateProperties/updateThermalTransportState.H"

if (kEffClosureModel_ == "Const")
{
    kappaEff  = kEffConst;
    pAlphaEff = kEffConst/Intr_CpF;
}
else if (kEffClosureModel_ == "ConstAnisotropic")
{
    // kappaEff remains the wall-normal (radial) conductivity used by the
    // porousMediaHeatTransferWall boundary condition.
    kappaEff  = kEffRadial;
    pAlphaEff = kEffRadial/Intr_CpF;
}
else if (kEffClosureModel_ == "Standard")
{
    // Standard volume-fraction-weighted conductivity from the manuscript.
    kappaEff  = por*Intr_kappaF + (1-por)*Intr_kappaS;
    pAlphaEff = kappaEff/Intr_CpF;
}
else if (kEffClosureModel_ == "Bunnell")
{
    kappaEff  = Intr_kappaF*(5+0.061*Re);
    pAlphaEff = kappaEff/Intr_CpF;
}
else if (kEffClosureModel_ == "Yagi")
{
    kappaEff  = Intr_kappaF*((1.0 - POR)*Intr_kappaS/(Intr_kappaF + YagiPhi*Intr_kappaS) + (YagiAlphaBeta*Re*Pr));
    pAlphaEff = kappaEff/Intr_CpF;
}
else if (kEffClosureModel_ == "Wakao")
{
    // ratios
    const volScalarField Rsk = Intr_kappaS/Intr_kappaF; // ks/kf
    const volScalarField Rfk = Intr_kappaF/Intr_kappaS; // kf/ks

    kappaEffStag = Intr_kappaF; // fallback

    if (stagnantKeffClosureModel_ == "Krupiczka")
    {
        // a = 0.280 - 0.757*log(por) - 0.057*log(ks/kf)
	volScalarField logPOR = Foam::log10(POR);
        volScalarField logR   = Foam::log10(Rsk);


	volScalarField a = 0.280 - 0.757*logPOR - 0.057*logR;

        // keff0/kf = (ks/kf)^a  => keff0 = kf * (ks/kf)^a
        kappaEffStag = Intr_kappaF * Foam::pow(Rsk, a);
    }
    else if (stagnantKeffClosureModel_ == "Kunii")
    {
        volScalarField denom = KuniiPhi + KuniiGamma*Rfk;
        // keff0/kf = por + (1-eps)/(phi + gamma*(kf/ks))
        denom = Foam::max(denom, dimensionedScalar("small", denom.dimensions(), SMALL));
        //
        volScalarField keff0_over_kf = POR + (1.0 - POR)/denom;

        kappaEffStag = Intr_kappaF * keff0_over_kf;

    }
    else if (stagnantKeffClosureModel_ == "ZS")
    {
        const volScalarField A   = Intr_kappaF/Intr_kappaS;
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

        kappaEffStag = Intr_kappaF * factor;
    }
    else
    {
        FatalIOErrorInFunction(porousMediaDict)
            << "Unknown stagnantKeffModel: " << stagnantKeffClosureModel_ << nl
            << "Valid: Krupiczka, Kunii, ZS"
            << exit(FatalIOError);
    }

        kappaEffMix = Cmix * Intr_kappaF * Re * Pr;

        kappaEff  = kappaEffStag + kappaEffMix;
        pAlphaEff = kappaEff/Intr_CpF;
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
        (kEffRadial/Intr_CpF)*symmTensor::I
      + ((kEffAxial - kEffRadial)/Intr_CpF)
       *sqr(kEffAxialDirection);
}
}
