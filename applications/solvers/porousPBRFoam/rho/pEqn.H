if ((!mesh.steady() && !pimple.simpleRho()) || pimple.consistent())
{
    Intr_rhoF = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(Intr_psiF*Intr_pF);

const volScalarField rAU("rAU", 1.0/UEqn.A());
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
    fvc::interpolate(Intr_rhoF)*fvc::flux(HbyA)
  + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(Intr_rhoF, Extr_UF, phi, rhoUf))
);

fvc::makeRelative(phiHbyA, Intr_rhoF, Extr_UF);
MRF.makeRelative(fvc::interpolate(Intr_rhoF), phiHbyA);

bool adjustMass = false;

// Update the pressure BCs to ensure flux consistency
constrainPressure(Intr_pF, Intr_rhoF, Extr_UF, phiHbyA, rhorAAtUf, MRF);

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

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(Intr_rhoF) + Intr_psiF*correction(fvm::ddt(Intr_pF))
      + fvc::div(phiHbyA) + fvm::div(phid, Intr_pF)
     ==
        fvOptions(Intr_psiF, Intr_pF, Intr_rhoF.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, Intr_pF));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    if (mesh.steady())
    {
        adjustMass = adjustPhi(phiHbyA, Extr_UF, Intr_pF);
    }

    if (pimple.consistent())
    {
        phiHbyA += (rhorAAtUf - rhorAUf)*fvc::snGrad(Intr_pF)*mesh.magSf();
        HbyA += (rAAtU - rAU)*fvc::grad(Intr_pF);
    }

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(Intr_rhoF) + Intr_psiF*correction(fvm::ddt(Intr_pF))
      + fvc::div(phiHbyA)
     ==
        fvOptions(Intr_psiF, Intr_pF, Intr_rhoF.name())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAAtUf, Intr_pF));

        pEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        pEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

if (mesh.steady())
{
    #include "incompressible/continuityErrs.H"
}
else
{
    const bool limitedp = pressureControl.limit(Intr_pF);

    // Thermodynamic density update
    thermo.correctRho(Intr_psiF*Intr_pF - psip0);

    if (limitedp)
    {
        Intr_rhoF = thermo.rho();
    }

    #include "porousMedia/equations/rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}

// Explicitly relax pressure for momentum corrector
Intr_pF.relax();

Extr_UF = HbyA - rAAtU*fvc::grad(Intr_pF);
Extr_UF.correctBoundaryConditions();
fvOptions.correct(Extr_UF);
K = 0.5*magSqr(Extr_UF);

if (mesh.steady())
{
    pressureControl.limit(Intr_pF);
}

// For steady compressible closed-volume cases adjust the pressure level
// to obey overall mass continuity
if (adjustMass && !thermo.incompressible())
{
    Intr_pF += (initialMass - fvc::domainIntegrate(thermo.rho()))
        /fvc::domainIntegrate(Intr_psiF);
    Intr_pF.correctBoundaryConditions();
}

if (mesh.steady() || pimple.simpleRho() || adjustMass)
{
    Intr_rhoF = thermo.rho();
}

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, Intr_rhoF, Extr_UF, phi);

if ((mesh.steady() || pimple.simpleRho()) && !pimple.transonic())
{
    Intr_rhoF.relax();
}

if (thermo.dpdt())
{
    Intr_dpdtF = fvc::ddt(Intr_pF);

    if (mesh.moving())
    {
        Intr_dpdtF -= fvc::div(fvc::meshPhi(Intr_rhoF, Extr_UF), Intr_pF);
    }
}
