    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "psi/createFields.H"
    #include "psi/createFieldsPorousMedia.H" //Nam

    if (solverMode_ == "PMSCHEM")
    {
	#include "psi/createFieldsCatalyst.H" //Nam
    }
    #include "createMRF.H"
    #include "createFvOptions.H"

    #include "psi/createFieldRefs.H"
    #include "psi/createCatalyticSource.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "psi/setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "porousMedia/equations/rhoEqn.H"

        while (pimple.loop())
        {           
            #include "psi/updateForPorousMedia.H" //Nam
            #include "porousMedia/equations/UEqn.H"

            // Surface chemistry contributes only when its model is present.
            #include "psi/calculateSurfaceReactions.H"

            // The same transport equations are used in PM and PMSCHEM.
            #include "porousMedia/equations/YEqn.H"

            if (!isothermal)
            {
                #include "porousMedia/equations/EEqn.H"
            }
            else
            {
                #include "porousMedia/equations/EEqn_iso.H"
            }

            #include "psi/setQdot.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "psi/pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        Intr_rhoF   = thermo.rho();
        Intr_muF    = thermo.mu();   //Nam
        Intr_CpF    = thermo.Cp();   //Nam
        Intr_kappaF = thermo.kappa(); //Nam
        Intr_heGasF = thermo.he(); //for cat

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
