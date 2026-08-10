{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        (POR+SMALL)*fvm::ddt(rho, he) + (1-POR-SMALL)*fvc::ddt(rhoS*CpS/Cp, he)
      + mvConvection->fvmDiv(phi, he)
      + (POR+SMALL)*
        (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(pAlphaEff, he)
     ==
        (POR+SMALL)*reaction->Qdot()
      + FS2V*surfaceReaction->QdotFromSCHEM()
      + (POR+SMALL)*fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;

    Info<< "min/max(p) = "
        << min(p).value() << ", " << max(p).value() << endl;
}
