tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    reaction->correct();
    volScalarField Yt(0.0*Intr_YF[0]);

    forAll(Intr_YF, i)
    {
        if (i != inertIndex && composition.active(i))
        {
            volScalarField& Intr_YiF = Intr_YF[i];

            tmp<fvScalarMatrix> tRcat
            (
                catalyticSpeciesSource(Intr_YiF)
            );

            fvScalarMatrix YiEqn
            (
                (POR+SMALL)*fvm::ddt(Intr_rhoF, Intr_YiF)
              + mvConvection->fvmDiv(phi, Intr_YiF)
             ==
                (POR+SMALL)*fvm::laplacian(Intr_rhoF*DiEff[i], Intr_YiF)
              + (POR+SMALL)*reaction->R(Intr_YiF)
              + tRcat()
              + (POR+SMALL)*fvOptions(Intr_rhoF, Intr_YiF)
            );

            YiEqn.relax();

            fvOptions.constrain(YiEqn);

            YiEqn.solve("Yi");

            fvOptions.correct(Intr_YiF);

            Intr_YiF.max(0.0);
            Yt += Intr_YiF;
        }
    }

    Intr_YF[inertIndex] = scalar(1) - Yt;
    Intr_YF[inertIndex].max(0.0);
}
