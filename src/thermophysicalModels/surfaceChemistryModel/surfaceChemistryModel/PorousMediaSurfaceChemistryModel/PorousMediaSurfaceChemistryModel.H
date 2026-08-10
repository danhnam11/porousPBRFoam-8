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

Class
    Foam::PorousMediaSurfaceChemistryModel

Description
    Extends base chemistry model by adding a thermo package, and ODE functions.
    Introduces chemistry equation system and evaluation of chemical source
    terms.

SourceFiles
    PorousMediaSurfaceChemistryModelI.H
    PorousMediaSurfaceChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef PorousMediaSurfaceChemistryModel_H
#define PorousMediaSurfaceChemistryModel_H

#include "BasicSurfaceChemistryModel.H"
#include "SurfaceReactionList.H"
#include "ReactionList.H"
#include "ODESystem.H"
#include "volFields.H"
#include "fvScalarMatrix.H" //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                      Class PorousMediaSurfaceChemistryModel Declaration
\*---------------------------------------------------------------------------*/

template
<
    class GasReactionThermo,
    class SolidReactionThermo,
    class GasThermoType,
    class SolidThermoType
>
class PorousMediaSurfaceChemistryModel
:
    public BasicSurfaceChemistryModel<GasReactionThermo, SolidReactionThermo>,
    public ODESystem
{
    // Private Member Functions

        //- Integrate reaction rate over the time-step
        //  using the selected ODE solver
        bool integrateSurfaceReactionRate_;

        bool usingMaxIntegrationTime_;

        //  Max integration time, reading from
        scalar maxIntegrationTime_;

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        template<class DeltaTType>
        scalar solvePatches(const DeltaTType& deltaT);

protected:

    //typedef GasThermoType gthermoType;
    typedef SolidThermoType thermoType;


    // Protected data

        //- Reference to the field of specie mass fractions
        //const PtrList<volScalarField>& Y_;
        PtrList<volScalarField>& Y_;

        //- Reference to the field of solid (site) specie mass fractions
        //const PtrList<volScalarField>& Ys_;
        PtrList<volScalarField>& Ys_;

        //- Thermodynamic data of the gas species
        const PtrList<GasThermoType>& gasSpecieThermos_;

        //- Thermodynamic data of the solid (site) species
        const PtrList<SolidThermoType>& specieThermos_;

        //- SurfaceReactions
        const ReactionList<GasThermoType> gasReactions_;

        //- SurfaceReactions
        const SurfaceReactionList<GasThermoType, SolidThermoType> reactions_;

        //- Number of species
        label nGasSpecie_;
        label nSolidSpecie_;
        label nSpecie_;

        //- Number of reactions
        label nReaction_;

        //- Temperature below which the reaction rates are assumed 0
        scalar Treact_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RRg_;
        PtrList<volScalarField::Internal> RRs_;
        volScalarField::Internal sumRRg_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary rate-of-change of concentration field
        mutable scalarField dcdt_;
        
        scalar FcatGeo_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Boundary> RRgwall_;
        PtrList<volScalarField::Boundary> RRswall_;
        volScalarField::Boundary sumRRgwall_;

        bool onlySurfaceCatalyst_; 

        List<word> catalyticSurfaces_;
        List<word> catalyticZones_;

        List<scalar> siteDensityList_;
        List<label> siteNumberList_;


    // Protected Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();



public:

    //- Runtime type information
    TypeName("porousMediaSCHEM");


    // Constructors

        //- Construct from thermo
        //PorousMediaSurfaceChemistryModel(const ReactionThermo& thermo);
        PorousMediaSurfaceChemistryModel
        (
            GasReactionThermo& thermo,
            SolidReactionThermo& catThermo
        );

        //- Disallow default bitwise copy construction
        PorousMediaSurfaceChemistryModel(const PorousMediaSurfaceChemistryModel&);


    //- Destructor
    virtual ~PorousMediaSurfaceChemistryModel();


    // Member Functions

        //- The reactions
        inline const PtrList<SurfaceReaction<GasThermoType, SolidThermoType>>& reactions() const;

        inline const PtrList<Reaction<GasThermoType>>& gasReactions() const;

        //- Thermodynamic data of the species
        inline const PtrList<SolidThermoType>& specieThermos() const;
        inline const PtrList<GasThermoType>& gasSpecieThermos() const;

        //- The number of species
        virtual inline label nSpecie() const;
        virtual inline label nGasSpecie() const;
        virtual inline label nSolidSpecie() const;

        //- The number of reactions
        virtual inline label nReaction() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar Treact() const;

        //- Temperature below which the reaction rates are assumed 0
        inline scalar& Treact();

        //- dc/dt = omega, rate of change in concentration, for each species
        virtual void omega
        (
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalarField& dcdti
        ) const;


        //- Return the reaction rate for iSurfaceReaction and the reference
        //  species and characteristic times
        virtual scalar omegaI
        (
            label iSurfaceReaction,
            const scalar p,
            const scalar T,
            const scalarField& c,
            const label li,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        //- Calculates the reaction rates
        virtual void calculate();

        // Chemistry model functions (overriding abstract functions in
        // basicCatalystModel.H)

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RRg
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RRg
            (
                const label i
            );

            //- Return const access to the chemical source terms for specie, i
            inline const volScalarField::Internal& RR
            (
                const label i
            ) const;

            //- Return non const access to chemical source terms [kg/m^3/s]
            virtual volScalarField::Internal& RR
            (
                const label i
            );


            //- Return reaction rate of the speciei in reactionI
            virtual tmp<volScalarField::Internal> calculateRR
            (
                const label reactionI,
                const label speciei
            ) const;

            //- Return reaction rate of reactionI
            virtual tmp<volScalarField::Internal> reactionRate
            (
                const label reactionI
            ) const;

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;


        // ODE functions (overriding abstract functions in ODE.H)

            //- Number of ODE's to solve
            inline virtual label nEqns() const;

            virtual void derivatives
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt
            ) const;

            virtual void jacobian
            (
                const scalar t,
                const scalarField& c,
                const label li,
                scalarField& dcdt,
                scalarSquareMatrix& J
            ) const;

            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& c,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const = 0;



        //- Calculates the surface reaction rates
        virtual void solveSurfaceReactions();

        //- Fuel consumption rate matrix from surf. reactions, 
	//  i.e. source term for fuel equation
        virtual tmp<fvScalarMatrix> RFromSCHEM(volScalarField& gasYi) const;

        //- Heat release rate from surface reactions [kg/m/s^3]
        virtual tmp<volScalarField> QdotFromSCHEM() const;

        //- Return species flux on patchi at catalytic surfaces
        virtual scalarField speciesFluxAtCatalyticSurfaces
        (
            const word speciei,
            const label patchi
        ) const;

        //- Return heat flux on patchi at catalytic surfaces
        virtual scalarField heatFluxAtCatalyticSurfaces
        (
            const label patchi
        ) const;

        //- Return momentum flux on patchi at catalytic surfaces
        virtual scalarField momentumFluxAtCatalyticSurfaces
        (
            const label patchi
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const PorousMediaSurfaceChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PorousMediaSurfaceChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PorousMediaSurfaceChemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
