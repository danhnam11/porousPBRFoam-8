// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    1./(POR+SMALL)*(fvm::ddt(rho, U) + 1./(POR+SMALL)*fvm::div(phi, U))
  + MRF.DDt(rho, U)
  + turbulence->divDevTau(U)
  + fvm::Sp(porousMomentumCoeff, U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == (-fvc::grad(p)));

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
