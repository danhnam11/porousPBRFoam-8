/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------


\*---------------------------------------------------------------------------*/

#include "speciesFluxFromCatalyticSurfacesFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicSurfaceChemistryModel.H"
#include "psiReactionThermophysicalTransportModel.H"
#include "basicThermo.H"
#include "psiReactionThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::
speciesFluxFromCatalyticSurfacesFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::
speciesFluxFromCatalyticSurfacesFvPatchScalarField
(
    const speciesFluxFromCatalyticSurfacesFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::
speciesFluxFromCatalyticSurfacesFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::
speciesFluxFromCatalyticSurfacesFvPatchScalarField
(
    const speciesFluxFromCatalyticSurfacesFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::
speciesFluxFromCatalyticSurfacesFvPatchScalarField
(
    const speciesFluxFromCatalyticSurfacesFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();
    const word speciei = IOobject::member(internalField().name());

    IOdictionary thermDict
    (
        IOobject
        (
            "thermophysicalProperties",
            this->db().time().constant(),
            this->db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    bool usingDetailedTransportModel
    (
        thermDict.lookupOrDefault("usingDetailedTransportModel", false)
    ); 

    bool soretEffect
    (
        thermDict.lookupOrDefault("soretEffect", false)
    ); 

    const psiReactionThermo& thermo = db().lookupObject<psiReactionThermo>
    (
        IOobject::groupName(basicThermo::dictName, IOobject::group(internalField().name()))
    );
    
    const basicSurfaceChemistryModel& surfaceReaction = db().lookupObject<basicSurfaceChemistryModel>
    (
        IOobject::groupName("surfaceChemistryProperties", IOobject::group(internalField().name()))
    );
   
    const psiReactionThermophysicalTransportModel& thermophysicalTransport = 
    db().lookupObject<psiReactionThermophysicalTransportModel>
    (
        IOobject::groupName("thermophysicalTransport", IOobject::group(internalField().name()))
    );

    // species flux from surface reactions, [kg/(m2.s)]
    const scalarField Rgw = 
        surfaceReaction.speciesFluxAtCatalyticSurfaces(speciei, patchi);            

    if (usingDetailedTransportModel)
    {
        label i = thermo.composition().species()[speciei];
        const wordList& speciesNames = thermo.composition().species();
        const volScalarField& Dimix = thermo.Dimix(i);        
        const scalarField rhow = thermo.rho(patchi);

        if (soretEffect)
        {

            const PtrList<volScalarField>& Y = thermo.composition().Y();
            const volScalarField T       = thermo.T();
            const volScalarField Wmix    = thermo.W();
            const volScalarField& Yi     = thermo.composition().Y(i);
            const volScalarField& DimixT = thermo.DimixT(i);
            
            const scalarField Tw         = T.boundaryField()[patchi];
            const scalarField TsnGrad    = T.boundaryField()[patchi].snGrad();
            const scalarField Wmixw      = Wmix.boundaryField()[patchi];
            const scalarField WmixsnGrad = Wmix.boundaryField()[patchi].snGrad();
            const scalarField Yw         = Yi.boundaryField()[patchi];
            const scalarField Dimixw     = Dimix.boundaryField()[patchi].patchInternalField();
            const scalarField DimixTw    = DimixT.boundaryField()[patchi].patchInternalField();
            
            // rho*D at wall-adjacent cell
            const scalarField rhoDimixw  = rhow*Dimixw;
            
            // Sum of species fluxes generated by surface reactions 
            // over all gas-phase species at this boundary
	    /*
            scalarField sumFlux(patch().size(), 0.0);
            forAll(Y, j)
            {
                sumFlux += surfaceReaction.speciesFluxAtCatalyticSurfaces
                (
                    speciesNames[j],
                    patchi
                );
            }
            */
            // Individual contributions
            const scalarField ReactionTerm   = Rgw/rhoDimixw;
            const scalarField WmixTerm       = Yw*WmixsnGrad/Wmixw;
            const scalarField SoretTerm      = DimixTw*TsnGrad/(rhoDimixw*Tw);
            //const scalarField ConvectiveTerm = Yw*sumFlux/(rhoDimixw);
            
            //gradient() = ReactionTerm - WmixTerm - SoretTerm - ConvectiveTerm; // [1/m]
            gradient() = ReactionTerm - WmixTerm - SoretTerm; // [1/m]
        }
        else
        {
            // without soret effect
            const fvPatchScalarField& DimixPf = Dimix.boundaryField()[patchi];
            const scalarField DimixBf(DimixPf);
            const scalarField rhoDimixBf = rhow*DimixBf;
            
            gradient() = Rgw/rhoDimixBf; // [1/m]
        }
    }
    else
    {
        // using unity Le and Sc assumptions
        const scalarField rhoDw = thermophysicalTransport.alphaEff(patchi);

        gradient() = Rgw/rhoDw; // [1/m]
    }

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::speciesFluxFromCatalyticSurfacesFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        speciesFluxFromCatalyticSurfacesFvPatchScalarField
    );
}

// ************************************************************************* //
