{
/*
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        (POR+SMALL)*fvm::ddt(rho, he) + (1-POR-SMALL)*fvc::ddt(rhoS*CpS/Cp, he)
      + mvConvection->fvmDiv(phi, he)
      + (POR+SMALL)*fvc::ddt(rho, K) + fvc::div(phi, K)  // for Porous media
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
      //+ thermophysicalTransport->divq(he) //original
      - fvm::laplacian(pAlphaEff, he) //ndn // for PORous media

     ==
        (POR+SMALL)*reaction->Qdot()
      //+ (POR+SMALL)*FS2V*catReaction->QdotFromSurf() //ndn for cat
      + FS2V*catReaction->QdotFromSurf() // F_{S2V}* dot{Q}^{cat}
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
*/
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
