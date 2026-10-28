// Isothermal operation: retain the prescribed temperature field.
    thermo.correct();

    Info<< "min/max(T) = "
        << min(Intr_TF).value() << ", " << max(Intr_TF).value() << endl;

    Info<< "min/max(p) = "
        << min(Intr_pF).value() << ", " << max(Intr_pF).value() << endl;
