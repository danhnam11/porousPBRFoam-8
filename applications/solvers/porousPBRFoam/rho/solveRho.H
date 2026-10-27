    #include "rho/checkSolverMode.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "rho/createFields.H"
    #include "rho/createFieldsPorousMedia.H" //Nam
    #include "rho/createFieldRefs.H"
    #include "rho/createCatalyticSource.H"
    #include "createRhoUfIfPresent.H"

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
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, Intr_rhoF, Extr_UF))
            );
        }

        if (LTS)
        {
            #include "rho/setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                // Store momentum to set rhoUf for introduced faces.
                autoPtr<volVectorField> rhoU;
                if (rhoUf.valid())
                {
                    rhoU = new volVectorField("rhoU", Intr_rhoF*Extr_UF);
                }

                // Do any mesh changes
                mesh.update();

                if (mesh.changing())
                {
                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & rhoUf();

                        #include "rho/correctPhi.H"

                        // Make the fluxes relative to the mesh-motion
                        fvc::makeRelative(phi, Intr_rhoF, Extr_UF);
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
            }

            if (pimple.firstPimpleIter() && !pimple.simpleRho())
            {
                #include "porousMedia/equations/rhoEqn.H"
            }

            #include "rho/updateForPorousMedia.H" //Nam

            #include "porousMedia/equations/UEqn.H"
            #include "porousMedia/equations/YEqn.H"
            if (!isothermal)
            {
                #include "porousMedia/equations/EEqn.H"
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "rho/pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        Intr_rhoF = thermo.rho();
        Intr_muF    = thermo.mu();   //Nam
        Intr_CpF    = thermo.Cp();   //Nam
        Intr_kappaF = thermo.kappa(); //Nam

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
