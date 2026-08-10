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
    Foam::porousMediaHeatTransferWallFvPatchScalarField

Description
    This boundary condition applies a heat flux condition to temperature
    on an external wall in one of three modes:

      - fixed power: supply Q
      - fixed heat flux: supply q
      - fixed heat transfer coefficient: supply h and Ta

    where:
    \vartable
        Q  | Power [W]
        q  | Heat flux [W/m^2]
        h  | Heat transfer coefficient [W/m^2/K]
        Ta | Ambient temperature [K]
    \endvartable

    In fixed-power and fixed-heat-flux modes, the conductive heat flux is
    evaluated using the porous-medium effective thermal conductivity,
    \c kappaEff.

    In fixed-heat-transfer-coefficient mode, the wall heat-transfer
    coefficient is evaluated using the closure selected through the
    \c WallHeatTransferClosureModel entry in \c porousMediaProperties.
    The available closures are \c Const, \c Yagi, \c Li, and \c Bey.

    Optional thin-layer thermal resistances can be specified through the
    \c thicknessLayers and \c kappaLayers entries. The ambient temperature,
    \c Ta, is specified as a Foam::Function1 of time and is spatially uniform.

Usage
    \table
    Property     | Description                 | Required | Default value
    mode         | 'power', 'flux' or 'coefficient' | yes |
    Q            | Power [W]                   | for mode 'power'     |
    q            | Heat flux [W/m^2]           | for mode 'flux'     |
    h            | Heat transfer coefficient [W/m^2/K] | for mode 'coefficent' |
    Ta           | Ambient temperature [K]     | for mode 'coefficient' |
    thicknessLayers | Layer thicknesses [m] | no |
    kappaLayers  | Layer thermal conductivities [W/m/K] | no |
    relaxation   | Relaxation for the wall temperature | no | 1
    emissivity   | Surface emissivity for radiative flux to ambient | no | 0
    qr           | Name of the radiative field | no | none
    qrRelaxation | Relaxation factor for radiative field | no | 1
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            porousMediaHeatTransferWall;

        mode            coefficient;

        Ta              constant 300.0;
        h               uniform 10.0;
        thicknessLayers (0.1 0.2 0.3 0.4);
        kappaLayers     (1 2 3 4);

        value           $internalField;
    }
    \endverbatim

See also
    Foam::temperatureCoupledBase
    Foam::mixedFvPatchScalarField

References
    Wall heat-transfer closures:

    [1] S. Yagi and D. Kunii, "Studies on heat transfer near wall
        surfaces in packed beds," AIChE J. 6 (1960) 97-104.

    [2] C.-H. Li and B. A. Finlayson, "Heat transfer in packed beds--a
        reevaluation," Chem. Eng. Sci. 32 (1977) 1055-1066.

    [3] O. Bey and G. Eigenberger, "Gas flow and heat transfer through
        catalyst filled tubes," Int. J. Therm. Sci. 40 (2001) 152-164.

    Packed-bed reactor application:

    [4] A. Takahashi and T. Fujitani, "Kinetic-model-based design of
        industrial reactor for catalytic hydrogen production via ammonia
        decomposition," Chem. Eng. Res. Des. 165 (2021) 333-340.

    by Danh Nam and Jae Hun Lee, CCER Lab, UNIST, Ulsan, Korea
    advisor, Prof. Chun Sang Yoo 

SourceFiles
    porousMediaHeatTransferWallFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef porousMediaHeatTransferWallFvPatchScalarField_H
#define porousMediaHeatTransferWallFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "temperatureCoupledBase.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class porousMediaHeatTransferWallFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class porousMediaHeatTransferWallFvPatchScalarField
:
    public mixedFvPatchScalarField,
    public temperatureCoupledBase
{
public:

    // Public data

        //- Operation mode enumeration
        enum operationMode
        {
            fixedPower,
            fixedHeatFlux,
            fixedHeatTransferCoeff
        };

        static const NamedEnum<operationMode, 3> operationModeNames;


private:

    // Private Data

        //- Operation mode
        operationMode mode_;

        //- Heat power [W]
        scalar Q_;

        //- Heat flux [W/m^2]
        scalarField q_;

        //- Heat transfer coefficient [W/m^2K]
        scalarField h_;

        //- Ambient temperature [K]
        autoPtr<Function1<scalar>> Ta_;

        //- Relaxation for the wall temperature (thermal inertia)
        scalar relaxation_;

        //- Optional surface emissivity for radiative transfer to ambient
        scalar emissivity_;

        //- Cache qr for relaxation
        scalarField qrPrevious_;

        //- Relaxation for qr
        scalar qrRelaxation_;

        //- Name of the radiative heat flux
        const word qrName_;

        //- Thickness of layers
        scalarList thicknessLayers_;

        //- Conductivity of layers
        scalarList kappaLayers_;

public:

    //- Runtime type information
    TypeName("porousMediaHeatTransferWall");


    // Constructors

        //- Construct from patch and internal field
        porousMediaHeatTransferWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        porousMediaHeatTransferWallFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  porousMediaHeatTransferWallFvPatchScalarField
        //  onto a new patch
        porousMediaHeatTransferWallFvPatchScalarField
        (
            const porousMediaHeatTransferWallFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        porousMediaHeatTransferWallFvPatchScalarField
        (
            const porousMediaHeatTransferWallFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new porousMediaHeatTransferWallFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        porousMediaHeatTransferWallFvPatchScalarField
        (
            const porousMediaHeatTransferWallFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new porousMediaHeatTransferWallFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Allow manipulation of the boundary values
            virtual bool fixesValue() const
            {
                return false;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
