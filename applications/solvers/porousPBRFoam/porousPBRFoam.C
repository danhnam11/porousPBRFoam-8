/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
//for cat
#include "psiReactionCatThermo.H"
#include "psiReactionCatThermo.H"
#include "SurfaceReaction.H"
#include "SurfaceReactionList.H"
#include "BasicSurfaceChemistryModel.H"
//end of for cat

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldsPorousMedia.H" //Nam

    if (solverMode_ == "PMSCHEM")
    {
	#include "createFieldsCatalyst.H" //Nam
    }
    #include "createMRF.H"
    #include "createFvOptions.H"

    #include "createFieldRefs.H"

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
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "PMrhoEqn.H"

        while (pimple.loop())
        {           
            #include "updateForPorousMedia.H" //Nam
            #include "PMUEqn.H"

    	    if (solverMode_ == "PMSCHEM")
	    {
	        #include "calculateSurfaceReactions.H" //Nam
	        #include "PMSCHEM_YEqn.H"

		if (!isothermal)
    		{
        	    #include "PMSCHEM_EEqn.H"
    		}
    		else
    		{
        	    #include "PMSCHEM_EEqn_iso.H"
    		}
	        //#include "PMSCHEM_EEqn.H"
        	#include "PMSCHEM_setQdot.H"
    	    }

	    else // solverMode_ == "PM"
	    {
	        #include "PM_YEqn.H"
		if (!isothermal)
		{
		    #include "PM_EEqn.H"
        	    #include "PM_setQdot.H"
		}
		else
		{
		    #include "PM_EEqn_iso.H"
		}
	        //#include "PM_EEqn.H"
        	//#include "PM_setQdot.H"
	    }
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        rho   = thermo.rho();
        mu    = thermo.mu();   //Nam
        Cp    = thermo.Cp();   //Nam
        kappa = thermo.kappa(); //Nam
        heGas = thermo.he(); //for cat

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
