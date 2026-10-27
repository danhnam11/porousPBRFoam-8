    forAll(DiEff, i)
    {
        DiEff[i] =
            sqrt(POR + SMALL)
           *(thermophysicalTransport->alphaEff()/Intr_rhoF)
          + (mag(Extr_UF)*DP)/(12.0*(POR + SMALL));
    }
