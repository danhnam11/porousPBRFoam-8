Intr_rhoF = thermo.rho();

const volScalarField rAU(1.0/UEqn.A());
const surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(Intr_rhoF*rAU));

tmp<volScalarField> rAtU
(
    pimple.consistent()
  ? volScalarField::New("rAtU", 1.0/(1.0/rAU - UEqn.H1()))
  : tmp<volScalarField>(nullptr)
);
tmp<surfaceScalarField> rhorAtUf
(
    pimple.consistent()
  ? surfaceScalarField::New("rhoRAtUf", fvc::interpolate(Intr_rhoF*rAtU()))
  : tmp<surfaceScalarField>(nullptr)
);

const volScalarField& rAAtU = pimple.consistent() ? rAtU() : rAU;
const surfaceScalarField& rhorAAtUf =
    pimple.consistent() ? rhorAtUf() : rhorAUf;

volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), Extr_UF, Intr_pF));

if (pimple.nCorrPiso() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(Intr_rhoF)*fvc::flux(HbyA)
      + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(Intr_rhoF, Extr_UF, phi))
    )
);

MRF.makeRelative(fvc::interpolate(Intr_rhoF), phiHbyA);

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(Intr_psiF)/fvc::interpolate(Intr_rhoF))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(Intr_psiF*Intr_pF)*phiHbyA/fvc::interpolate(Intr_rhoF);

    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(Intr_pF)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(Intr_pF);
    }

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(Intr_psiF, Intr_pF)
          + fvc::div(phiHbyA) + fvm::div(phid, Intr_pF)
          - fvm::laplacian(rhorAAtUf, Intr_pF)
         ==
            fvOptions(Intr_psiF, Intr_pF, Intr_rhoF.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi == phiHbyA + pEqn.flux();
        }
    }
}
else
{
    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(Intr_pF)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(Intr_pF);
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(Intr_pF, Intr_rhoF, Extr_UF, phiHbyA, rhorAAtUf, MRF);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(Intr_psiF, Intr_pF)
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAAtUf, Intr_pF)
         ==
            fvOptions(Intr_psiF, Intr_pF, Intr_rhoF.name())
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "porousMedia/equations/rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
Intr_pF.relax();

Extr_UF = HbyA - rAAtU*fvc::grad(Intr_pF);
Extr_UF.correctBoundaryConditions();
fvOptions.correct(Extr_UF);
K = 0.5*magSqr(Extr_UF);

if (pressureControl.limit(Intr_pF))
{
    Intr_pF.correctBoundaryConditions();
    Intr_rhoF = thermo.rho();
}

if (thermo.dpdt())
{
    Intr_dpdtF = fvc::ddt(Intr_pF);
}
