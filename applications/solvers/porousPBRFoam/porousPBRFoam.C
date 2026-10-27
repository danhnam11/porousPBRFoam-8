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
    porousPBRFoam_integrate

Description
    Unified porous packed-bed solver. The thermodynamics entry selects
    the psi or rho numerical path at run time.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "psi/psiTransportCompatibility.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "psiReactionCatThermo.H"
#include "SurfaceReaction.H"
#include "SurfaceReactionList.H"
#include "BasicSurfaceChemistryModel.H"
#include "dynamicFvMesh.H"
#include "rhoReactionThermo.H"
#include "rhoReactionThermophysicalTransportModel.H"
#include "CorrectPhi.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Porous packed-bed solver with hePsiThermo/heRhoThermo selection."
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "initContinuityErrs.H"

    #include "readThermoType.H"

    if (thermoType == "hePsiThermo")
    {
        Info<< "Selecting psi thermodynamics for porousPBRFoam_integrate" << nl;
        #include "psi/solvePsi.H"
    }
    else if (thermoType == "heRhoThermo")
    {
        Info<< "Selecting rho thermodynamics for porousPBRFoam_integrate" << nl;
        #include "rho/solveRho.H"
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported thermoType.type: " << thermoType << nl
            << "Expected hePsiThermo or heRhoThermo."
            << exit(FatalError);
    }
    return 0;
}

// ************************************************************************* //
