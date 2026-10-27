{
    volScalarField& Intr_heF = thermo.he();

    fvScalarMatrix EEqn
    (
        (POR+SMALL)*fvm::ddt(Intr_rhoF, Intr_heF) + (1-POR-SMALL)*fvc::ddt(Intr_rhoS*Intr_CpS/Intr_CpF, Intr_heF)
      + mvConvection->fvmDiv(phi, Intr_heF)
      + (POR+SMALL)*
        (
            Intr_heF.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(Intr_rhoF), Extr_UF),
                Intr_pF,
                "div(phiv,p)"
            )
          : -Intr_dpdtF
        )
      - fvm::laplacian(pAlphaEffTensor, Intr_heF)
     ==
        (POR+SMALL)*reaction->Qdot()
      + QdotSurface
      + (POR+SMALL)*fvOptions(Intr_rhoF, Intr_heF)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(Intr_heF);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(Intr_TF).value() << ", " << max(Intr_TF).value() << endl;

    Info<< "min/max(p) = "
        << min(Intr_pF).value() << ", " << max(Intr_pF).value() << endl;
}
