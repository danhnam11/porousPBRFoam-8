// Solve the Momentum equation

MRF.correctBoundaryVelocity(Extr_UF);

tmp<fvVectorMatrix> tUEqn
(
    1./(POR+SMALL)*(fvm::ddt(Intr_rhoF, Extr_UF) + 1./(POR+SMALL)*fvm::div(phi, Extr_UF))
  + MRF.DDt(Intr_rhoF, Extr_UF)
  + turbulence->divDevTau(Extr_UF)
  + fvm::Sp(porousMomentumCoeff, Extr_UF)
 ==
    fvOptions(Intr_rhoF, Extr_UF)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == (-fvc::grad(Intr_pF)));

    fvOptions.correct(Extr_UF);
    K = 0.5*magSqr(Extr_UF);
}
