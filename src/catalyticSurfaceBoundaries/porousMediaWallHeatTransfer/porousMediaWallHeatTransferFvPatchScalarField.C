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

    This boundary condition has been extended for heat transfer between a
    reactor wall and a porous packed bed. The implementation:
        1. uses the porous-medium effective thermal conductivity, kappaEff,
           to evaluate the conductive heat flux at the wall; and
        2. evaluates the wall heat-transfer coefficient, h_, using the
           closure selected in porousMediaProperties.

    The available wall heat-transfer closures are Const, Yagi, Li, and Bey.

    References for the wall heat-transfer closures:
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
   
\*---------------------------------------------------------------------------*/

#include "porousMediaWallHeatTransferFvPatchScalarField.H"
#include "volFields.H"
#include "physicoChemicalConstants.H"
#include "addToRunTimeSelectionTable.H"
#include "thermophysicalTransportModel.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include <sstream>
#include <iomanip>
#include <vector>
#include "DynamicList.H"
#include <cmath>

// Wall-condition input validation stays local to this library. Solver closure
// validation is kept separately under applications/solvers/porousPBRFoam.
namespace Foam
{
namespace
{
inline void requirePorousWallInputs
(
    const dictionary& dict,
    const string& model,
    const char* names,
    const string& description
)
{
    IStringStream requiredStream(names);
    const wordList required(requiredStream);
    DynamicList<word> missing;
    forAll(required, i)
    {
        if (!dict.found(required[i])) missing.append(required[i]);
    }
    if (missing.size())
    {
        // Use fixed-width, left-aligned columns in a standard text stream.
        // Determine widths from every row before emitting the table.
        std::vector<std::string> symbols(required.size(), "-");
        std::vector<std::string> units(required.size(), "-");
        forAll(required, i)
        {
            const std::string keyword(required[i].c_str());
            std::istringstream parameters(description.c_str());
            std::string parameter;
            while (std::getline(parameters, parameter, ';'))
            {
                const std::string::size_type first =
                    parameter.find_first_not_of(' ');
                if (first == std::string::npos) continue;
                parameter.erase(0, first);
                const std::string prefix = keyword + " = ";
                if (parameter.compare(0, prefix.size(), prefix) != 0) continue;
                const std::string body = parameter.substr(prefix.size());
                const std::string::size_type open = body.find('[');
                const std::string::size_type close = body.find(']', open);
                if (open != std::string::npos && close != std::string::npos)
                {
                    symbols[i] = body.substr(0, open);
                    const std::string::size_type last =
                        symbols[i].find_last_not_of(' ');
                    symbols[i].erase(last == std::string::npos ? 0 : last + 1);
                    units[i] = body.substr(open + 1, close - open - 1);
                }
            }
            if (keyword == "AxialDirection")
                symbols[i] = "e_z";
            else if (keyword == "StagnantThermalConductivityClosureModel")
                symbols[i] = "k_eff^0 model";
            else if (keyword == "ThermalConductivityClosureModel")
                symbols[i] = "k_eff model";
        }

        std::size_t keywordWidth = 28;
        std::size_t symbolWidth = 16;
        std::size_t unitWidth = 10;
        forAll(required, i)
        {
            if (required[i].size() > keywordWidth)
                keywordWidth = required[i].size();
            if (symbols[i].size() > symbolWidth) symbolWidth = symbols[i].size();
            if (units[i].size() > unitWidth) unitWidth = units[i].size();
        }
        const std::string selection(model.c_str());
        const std::string::size_type separator = selection.find(" = ");
        const std::string modelName = separator == std::string::npos
            ? selection : selection.substr(separator + 3);
        const std::string entryName = separator == std::string::npos
            ? "-" : selection.substr(0, separator);

        std::ostringstream message;
        message << "\nMissing required closure-model inputs\n\n"
            << "Model       : " << modelName << '\n'
            << "Selection   : " << entryName << '\n'
            << "Dictionary  : constant/porousMediaProperties\n\n"
            << std::left
            << std::setw(keywordWidth) << "Keyword" << "  "
            << std::setw(symbolWidth) << "Symbol" << "  "
            << std::setw(unitWidth) << "Unit" << "  "
            << "Status\n"
            << std::string(keywordWidth + symbolWidth + unitWidth + 14, '-')
            << '\n';
        forAll(required, i)
        {
            message << std::setw(keywordWidth) << required[i].c_str() << "  "
                << std::setw(symbolWidth) << symbols[i] << "  "
                << std::setw(unitWidth) << units[i] << "  "
                << (dict.found(required[i]) ? "provided" : "MISSING") << '\n';
        }
        message
            << std::string(keywordWidth + symbolWidth + unitWidth + 14, '-')
            << "\n\nAction      : Specify the entries marked MISSING.\n";

        FatalIOErrorInFunction(dict)
            << message.str().c_str() << exit(FatalIOError);
    }
}

inline word porousWallModelSelection
(
    const dictionary& dict,
    const word& keyword,
    const char* names
)
{
    IStringStream validStream(names);
    const wordList valid(validStream);
    if (!dict.found(keyword))
    {
        OStringStream message;
        message << nl
            << "Porous-media closure: missing model selection" << nl
            << "--------------------------------------------" << nl
            << "Entry      : " << keyword << nl
            << "Dictionary : constant/porousMediaProperties" << nl
            << "Valid options: ";
        forAll(valid, i)
        {
            if (i) message << ", ";
            message << valid[i];
        }
        message << nl << "Action: specify the model-selection entry." << nl;
        FatalIOErrorInFunction(dict)
            << message.str().c_str() << exit(FatalIOError);
    }
    const word selected(dict.lookup(keyword));
    bool found = false;
    forAll(valid, i) found = found || selected == valid[i];
    if (!found)
    {
        OStringStream message;
        message << nl
            << "Porous-media closure: unknown model" << nl
            << "----------------------------------" << nl
            << "Entry      : " << keyword << nl
            << "Value      : " << selected << nl
            << "Dictionary : constant/porousMediaProperties" << nl
            << "Valid options: ";
        forAll(valid, i)
        {
            if (i) message << ", ";
            message << valid[i];
        }
        message << nl << "Action: select one of the models listed above." << nl;
        FatalIOErrorInFunction(dict)
            << message.str().c_str() << exit(FatalIOError);
    }
    return selected;
}

inline void checkPorousWallScalarInput
(
    const dictionary& dict,
    const word& keyword,
    const string& model,
    const string& description,
    const bool allowZero = false
)
{
    const scalar value = readScalar(dict.lookup(keyword));
    if (!std::isfinite(value) || (allowZero ? value < 0 : value <= 0))
    {
        FatalIOErrorInFunction(dict)
            << nl << "Porous-media closure: invalid value" << nl
            << "----------------------------------" << nl
            << "Selection  : " << model.c_str() << nl
            << "Entry      : " << keyword << nl
            << "Value      : " << value << nl
            << "Parameter  : " << description.c_str() << nl
            << "Required   : finite and "
            << (allowZero ? "non-negative." : "strictly positive.") << nl
            << "Dictionary : constant/porousMediaProperties" << nl
            << "Action: correct the value of " << keyword << '.' << nl
            << exit(FatalIOError);
    }
}

inline void validatePorousWallInputs(const dictionary& dict)
{
    const word model = porousWallModelSelection
    (dict, "WallHeatTransferClosureModel", "(Const Yagi Li Bey)");
    const string context("WallHeatTransferClosureModel = " + model);
    const char* names = "(SolidThermalConductivity SuperficialVelocity ParticleDiameter ThermalConductivityClosureModel)";
    const char* description = "SolidThermalConductivity = k_s [W/m/K]; "
        "SuperficialVelocity = |<u>| [m/s]; ParticleDiameter = d_p [m]; "
        "ThermalConductivityClosureModel selects k_eff";
    if (model == "Const")
    {
        names = "(SolidThermalConductivity heatTransferCoefficient ThermalConductivityClosureModel)";
        description = "SolidThermalConductivity = k_s [W/m/K]; "
            "heatTransferCoefficient = h_w [W/m2/K]; "
            "ThermalConductivityClosureModel selects k_eff";
    }
    else if (model == "Bey")
    {
        names = "(SolidThermalConductivity SuperficialVelocity ParticleDiameter ReactorDiameter ThermalConductivityClosureModel)";
        description = "SolidThermalConductivity = k_s [W/m/K]; "
            "SuperficialVelocity = |<u>| [m/s]; ParticleDiameter = d_p [m]; "
            "ReactorDiameter = d_t [m]; ThermalConductivityClosureModel = Wakao "
            "provides k_eff^0";
    }
    requirePorousWallInputs(dict, context, names, description);
    checkPorousWallScalarInput(dict, "SolidThermalConductivity", context, "k_s [W/m/K]");
    if (model == "Const")
        checkPorousWallScalarInput(dict, "heatTransferCoefficient", context, "h_w [W/m2/K]", true);
    else
    {
        checkPorousWallScalarInput(dict, "SuperficialVelocity", context, "|<u>| [m/s]", true);
        checkPorousWallScalarInput(dict, "ParticleDiameter", context, "d_p [m]");
    }
    if (model == "Bey")
    {
        checkPorousWallScalarInput(dict, "ReactorDiameter", context, "d_t [m]");
        const word thermal(dict.lookup("ThermalConductivityClosureModel"));
        if (thermal != "Wakao")
        {
            FatalIOErrorInFunction(dict)
                << context << " requires ThermalConductivityClosureModel = Wakao "
                << "to provide the stagnant effective conductivity k_eff^0 "
                << "(kappaEffStag). Selected: " << thermal
                << exit(FatalIOError);
        }
        if (readScalar(dict.lookup("ReactorDiameter"))
         <= readScalar(dict.lookup("ParticleDiameter")))
        {
            FatalIOErrorInFunction(dict)
                << context << " requires ReactorDiameter (d_t) > "
                << "ParticleDiameter (d_p)." << exit(FatalIOError);
        }
    }
    porousWallModelSelection
    (
        dict, "ThermalConductivityClosureModel",
        "(Const ConstAnisotropic Standard Bunnell Yagi Wakao)"
    );
}

} // End anonymous namespace
} // End namespace Foam

using Foam::constant::physicoChemical::sigma;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        porousMediaWallHeatTransferFvPatchScalarField::operationMode,
        3
    >::names[] =
    {
        "power",
        "flux",
        "coefficient"
    };
}

const Foam::NamedEnum
<
    Foam::porousMediaWallHeatTransferFvPatchScalarField::operationMode,
    3
> Foam::porousMediaWallHeatTransferFvPatchScalarField::operationModeNames;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousMediaWallHeatTransferFvPatchScalarField::
porousMediaWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch()),
    mode_(fixedHeatFlux),
    Q_(0),
    Ta_(),
    relaxation_(1),
    emissivity_(0),
    qrRelaxation_(1),
    qrName_("undefined-qr"),
    thicknessLayers_(),
    kappaLayers_()
{
    refValue() = 0;
    refGrad() = 0;
    valueFraction() = 1;
}


Foam::porousMediaWallHeatTransferFvPatchScalarField::
porousMediaWallHeatTransferFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    mode_(operationModeNames.read(dict.lookup("mode"))),
    Q_(0),
    Ta_(),
    relaxation_(dict.lookupOrDefault<scalar>("relaxation", 1)),
    emissivity_(dict.lookupOrDefault<scalar>("emissivity", 0)),
    qrRelaxation_(dict.lookupOrDefault<scalar>("qrRelaxation", 1)),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(),
    kappaLayers_()
{
    switch (mode_)
    {
        case fixedPower:
        {
            dict.lookup("Q") >> Q_;

            break;
        }
        case fixedHeatFlux:
        {
            q_ = scalarField("q", dict, p.size());

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_ = scalarField("h", dict, p.size());
            Ta_ = Function1<scalar>::New("Ta", dict);

            if (dict.found("thicknessLayers"))
            {
                dict.lookup("thicknessLayers") >> thicknessLayers_;
                dict.lookup("kappaLayers") >> kappaLayers_;
            }

            break;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (qrName_ != "none")
    {
        if (dict.found("qrPrevious"))
        {
            qrPrevious_ = scalarField("qrPrevious", dict, p.size());
        }
        else
        {
            qrPrevious_.setSize(p.size(), 0);
        }
    }

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0;
        valueFraction() = 1;
    }
}


Foam::porousMediaWallHeatTransferFvPatchScalarField::
porousMediaWallHeatTransferFvPatchScalarField
(
    const porousMediaWallHeatTransferFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    temperatureCoupledBase(patch(), ptf),
    mode_(ptf.mode_),
    Q_(ptf.Q_),
    Ta_(ptf.Ta_, false),
    relaxation_(ptf.relaxation_),
    emissivity_(ptf.emissivity_),
    qrRelaxation_(ptf.qrRelaxation_),
    qrName_(ptf.qrName_),
    thicknessLayers_(ptf.thicknessLayers_),
    kappaLayers_(ptf.kappaLayers_)
{
    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            mapper(q_, ptf.q_);
            break;
        }
        case fixedHeatTransferCoeff:
        {
            mapper(h_, ptf.h_);
            break;
        }
    }

    if (qrName_ != "none")
    {
        mapper(qrPrevious_, ptf.qrPrevious_);
    }
}


Foam::porousMediaWallHeatTransferFvPatchScalarField::
porousMediaWallHeatTransferFvPatchScalarField
(
    const porousMediaWallHeatTransferFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf),
    temperatureCoupledBase(tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


Foam::porousMediaWallHeatTransferFvPatchScalarField::
porousMediaWallHeatTransferFvPatchScalarField
(
    const porousMediaWallHeatTransferFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF),
    temperatureCoupledBase(patch(), tppsf),
    mode_(tppsf.mode_),
    Q_(tppsf.Q_),
    q_(tppsf.q_),
    h_(tppsf.h_),
    Ta_(tppsf.Ta_, false),
    relaxation_(tppsf.relaxation_),
    emissivity_(tppsf.emissivity_),
    qrPrevious_(tppsf.qrPrevious_),
    qrRelaxation_(tppsf.qrRelaxation_),
    qrName_(tppsf.qrName_),
    thicknessLayers_(tppsf.thicknessLayers_),
    kappaLayers_(tppsf.kappaLayers_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousMediaWallHeatTransferFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            m(q_, q_);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            m(h_, h_);

            break;
        }
    }

    if (qrName_ != "none")
    {
        m(qrPrevious_, qrPrevious_);
    }
}


void Foam::porousMediaWallHeatTransferFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const porousMediaWallHeatTransferFvPatchScalarField& tiptf =
        refCast<const porousMediaWallHeatTransferFvPatchScalarField>(ptf);

    switch (mode_)
    {
        case fixedPower:
        {
            break;
        }
        case fixedHeatFlux:
        {
            q_.rmap(tiptf.q_, addr);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            h_.rmap(tiptf.h_, addr);

            break;
        }
    }

    if (qrName_ != "none")
    {
        qrPrevious_.rmap(tiptf.qrPrevious_, addr);
    }
}


void Foam::porousMediaWallHeatTransferFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Tp(*this);

    // Store current valueFraction and refValue for relaxation
    const scalarField valueFraction0(valueFraction());
    const scalarField refValue0(refValue());

    scalarField qr(Tp.size(), 0);
    if (qrName_ != "none")
    {
        qr =
            qrRelaxation_
           *patch().lookupPatchField<volScalarField, scalar>(qrName_)
          + (1 - qrRelaxation_)*qrPrevious_;

        qrPrevious_ = qr;
    }

    switch (mode_)
    {
        case fixedPower:
        {
            // Read dictionary file named "porousMediaProperties"
            // placed at constant directory as input

            IOdictionary porousMediaDict
            (
                IOobject
                (
                    "porousMediaProperties",
                    this->db().time().constant(),
                    this->db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const fvMesh& mesh_ = patch().boundaryMesh().mesh();
            const label patchi = patch().index();

            if (!mesh_.foundObject<volScalarField>("kappaEff"))
            {
                FatalErrorInFunction
                    << "Required volScalarField 'kappaEff' not found in mesh database.\n"
                    << "This boundary condition assumes that 'kappaEff' is computed\n"
                    << "and registered by the solver or a model before BC evaluation.\n"
                    << "Patch      : " << patch().name() << nl
                    << "Field      : " << this->internalField().name() << nl
                    << "Region     : " << mesh_.name() << nl
                    << exit(FatalError);
            }

            const volScalarField& kappaEff =
                mesh_.lookupObject<volScalarField>("kappaEff");

            const scalarField& kappaEffw =
                kappaEff.boundaryField()[patchi];

            // Use the effective thermal conductivity at the wall
            refGrad() = (Q_/gSum(patch().magSf()) + qr)/kappaEffw;
            refValue() = Tp;
            valueFraction() = 0;

            break;
        }
        case fixedHeatFlux:
        {
            // Read dictionary file named "porousMediaProperties"
            // placed at constant directory as input
            
            IOdictionary porousMediaDict
            (
                IOobject
                (
                    "porousMediaProperties",
                    this->db().time().constant(),
                    this->db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            const fvMesh& mesh_ = patch().boundaryMesh().mesh();
            const label patchi = patch().index();

	    if (!mesh_.foundObject<volScalarField>("kappaEff"))
 	    {
    		FatalErrorInFunction
        	    << "Required volScalarField 'kappaEff' not found in mesh database.\n"
        	    << "This boundary condition assumes that 'kappaEff' is computed\n"
        	    << "and registered by the solver or a model before BC evaluation.\n"
        	    << "Patch      : " << patch().name() << nl
        	    << "Field      : " << this->internalField().name() << nl
        	    << "Region     : " << mesh_.name() << nl
        	    << exit(FatalError);
	    }

            const volScalarField& kappaEff =
                mesh_.lookupObject<volScalarField>("kappaEff");

	    const scalarField& kappaEffw =
	        kappaEff.boundaryField()[patchi];

	    // Use the effective thermal conductivity at the wall
            refGrad() = (q_ + qr)/kappaEffw;
            refValue() = Tp;
            valueFraction() = 0;

            break;
        }
        case fixedHeatTransferCoeff:
        {
            scalar totalSolidRes = 0;
            if (thicknessLayers_.size())
            {
                forAll(thicknessLayers_, iLayer)
                {
                    const scalar l = thicknessLayers_[iLayer];
                    if (kappaLayers_[iLayer] > 0)
                    {
                        totalSolidRes += l/kappaLayers_[iLayer];
                    }
                }
            }

            // Read dictionary file named "porousMediaProperties"
            // placed at constant directory as input
            IOdictionary porousMediaDict
            (
                IOobject
                (
                    "porousMediaProperties",
                    this->db().time().constant(),
                    this->db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            // Validate model requirements before any per-entry lookup.
            validatePorousWallInputs(porousMediaDict);

            // Solid thermal conductivity [W/m-K]
            scalar ks(readScalar(porousMediaDict.lookup("SolidThermalConductivity")));                        
            
            // Select method to evaluate the heat transfer coefficient at walls
            word hwClosureModel_ = word(porousMediaDict.lookup("WallHeatTransferClosureModel"));

	    if
	    (
	    	hwClosureModel_ != "Const"
	    	&& hwClosureModel_ != "Yagi"
	        && hwClosureModel_ != "Li"
	        && hwClosureModel_ != "Bey"
	    )
            {
                FatalErrorInFunction
                    << "Unknown WallHeatTransferClosureModel \""
                    << hwClosureModel_
                    << "\" specified in porousMediaProperties." << nl
                    << "Valid options are:" << nl
                    << "  Const" << nl
                    << "  Yagi" << nl
                    << "  Li" << nl
                    << "  Bey" << nl
                    << nl
                    << "Please check the dictionary entry:" << nl
                    << "  WallHeatTransferClosureModel <modelName>;" << nl
                    << exit(FatalError);
            }

            bool useConstHW   = (hwClosureModel_   == "Const");
            bool useYagiHW    = (hwClosureModel_   == "Yagi");
            bool useLiHW      = (hwClosureModel_   == "Li");
            bool useBeyHW     = (hwClosureModel_   == "Bey");

	    word kEffClosureModel(porousMediaDict.lookup("ThermalConductivityClosureModel"));

            // Access the boundary patch and its adjacent cells
            const fvMesh& mesh_ = patch().boundaryMesh().mesh();
            const label patchi = patch().index();
            const Foam::fvPatch& currPatch = mesh_.boundary()[patchi];
            const labelUList& faceCells = currPatch.faceCells();
            scalar faceCellNum = faceCells.size();

            if (!mesh_.foundObject<volScalarField>("kappaEff"))
            {
                FatalErrorInFunction
                    << "Required volScalarField 'kappaEff' not found in mesh database.\n"
                    << "This boundary condition assumes that 'kappaEff' is computed\n"
                    << "and registered by the solver or a model before BC evaluation.\n"
                    << "Patch      : " << patch().name() << nl
                    << "Field      : " << this->internalField().name() << nl
                    << "Region     : " << mesh_.name() << nl
                    << exit(FatalError);
            }

            const volScalarField &kappaEff = mesh_.lookupObject<volScalarField>("kappaEff");

            scalarField kappaEff0w(patch().size(), 0.0);
        
            if (kEffClosureModel == "Wakao")
            {
                if (!mesh_.foundObject<volScalarField>("kappaEffStag"))
                {
                    FatalErrorInFunction
                        << "kEffClosureModel is Wakao, but 'kappaEffStag' was not found "
                        << "in the mesh database." << nl
                        << "Patch  : " << patch().name() << nl
                        << "Field  : " << this->internalField().name() << nl
                        << "Region : " << mesh_.name() << nl
                        << exit(FatalError);
                }
        
                const volScalarField& kappaEff0 =
                    mesh_.lookupObject<volScalarField>("kappaEffStag");
        
                kappaEff0w = kappaEff0.boundaryField()[patchi];
            }


            const scalarField& kappaEffw =
                kappaEff.boundaryField()[patchi];

            // Expand the scalar solid conductivity over the patch faces
            scalarField kappasw(faceCellNum, ks);

            // Access the local thermophysical and transport properties
            const thermophysicalTransportModel& ttm =
                db().lookupObject<thermophysicalTransportModel>
                (   
                    IOobject::groupName
                    (   
                        thermophysicalTransportModel::typeName,
                        internalField().group()
                    )
                );
            
            const compressibleMomentumTransportModel& turbModel =
                ttm.momentumTransport();
            
            const scalarField kappaw(ttm.kappaEff(patchi));
            
            const tmp<scalarField> tmuw = turbModel.mu(patchi);
            const scalarField& muw = tmuw();
            const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
            const scalarField& Tw = ttm.thermo().T().boundaryField()[patchi];
            const scalarField Cpw(ttm.thermo().Cp(Tw, patchi));

            // calculate Prandtl number
            const scalarField Pr(muw*Cpw/kappaw);


	    // Initialize geometric and operating parameters used by closures
	    scalar usup = 0.0;
	    scalar dp   = 0.0;
	    scalar dt   = 0.0;

            // update heat transfer coefficient as in eq. 31 as presented in Yagi et al.*
            // otherwise it is const number read from dictionary as default
	    if (useYagiHW || useLiHW)
	    {
		usup = readScalar(porousMediaDict.lookup("SuperficialVelocity"));
		dp   = readScalar(porousMediaDict.lookup("ParticleDiameter"));
	    }

            if (useBeyHW)
            {
                usup = readScalar(porousMediaDict.lookup("SuperficialVelocity"));
                dp   = readScalar(porousMediaDict.lookup("ParticleDiameter"));
                dt   = readScalar(porousMediaDict.lookup("ReactorDiameter"));
            }

	    // use Const wall heat transfer correlation
            if (useConstHW)
            {
                scalar hconst = readScalar(porousMediaDict.lookup("heatTransferCoefficient"));

                forAll(h_, faceI)
                {
                    h_[faceI] = hconst;
                }
            }
	    // use Yagi's wall heat transfer correlation
	    else if (useYagiHW)
	    {
	  	forAll(h_, faceI)
	  	{
	  	    scalar Ref = rhow[faceI]*dp*usup/muw[faceI];
                    scalar Prf = Pr[faceI];
	  	    h_[faceI]  = (kappaw[faceI]/dp)*(3+0.054*Ref*Prf);
	        }
	    }
	    // use Li & Finlayson's wall heat transfer correlation
            else if (useLiHW)
            {
                forAll(h_, faceI)
                {
                    scalar Ref = rhow[faceI]*dp*usup/muw[faceI];
                    h_[faceI]  = (kappaw[faceI]/dp)*(0.17*std::pow(Ref,0.79));
                }
            }
	    // use Bey & Eigenberger's wall heat transfer correlation
            else if (useBeyHW)
            {
                forAll(h_, faceI)
                {
                    scalar Ref = rhow[faceI]*dp*usup/muw[faceI];
                    scalar Prf = Pr[faceI];
                    h_[faceI]  = (kappaw[faceI]/dp)*
	            	(2.4*kappaEff0w[faceI]/kappaw[faceI]
		      + 0.054*(1-(dp/dt))*Ref*pow(Prf,(1/3)));
                }
            }
	    
            scalarField kappaDeltaCoeffs(faceCellNum, 1.0);

            // update effective thermal conductivity at the wall as in eq. 25
            // as presented in Yagi et al.* 
            // otherwise it is calculated by standard form (i.e., same as ANSYS)
            
            kappaDeltaCoeffs = kappaEffw*patch().deltaCoeffs();

            const scalar Ta = Ta_->value(this->db().time().timeOutputValue());
            const scalarField hp
            (
                1
               /(
                    1
                   /(
                        (emissivity_ > 0)
                      ? (
                            h_
                          + emissivity_*sigma.value()
                           *((pow3(Ta) + pow3(Tp)) + Ta*Tp*(Ta + Tp))
                        )()
                      : h_
                    ) + totalSolidRes
                )
            );

            const scalarField hpTa(hp*Ta);

            refGrad() = 0;

            forAll(Tp, i)
            {
                if (qr[i] < 0)
                {
                    const scalar hpmqr = hp[i] - qr[i]/Tp[i];

                    refValue()[i] = hpTa[i]/hpmqr;
                    valueFraction()[i] = hpmqr/(hpmqr + kappaDeltaCoeffs[i]);
                }
                else
                {
                    refValue()[i] = (hpTa[i] + qr[i])/hp[i];
                    valueFraction()[i] = hp[i]/(hp[i] + kappaDeltaCoeffs[i]);
                }
            }

            break;
        }
    }

    valueFraction() =
        relaxation_*valueFraction()
      + (1 - relaxation_)*valueFraction0;

    refValue() = relaxation_*refValue() + (1 - relaxation_)*refValue0;

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        const scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::porousMediaWallHeatTransferFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);

    writeEntry(os, "mode", operationModeNames[mode_]);
    temperatureCoupledBase::write(os);

    switch (mode_)
    {
        case fixedPower:
        {
            writeEntry(os, "Q", Q_);

            break;
        }
        case fixedHeatFlux:
        {
            writeEntry(os, "q", q_);

            break;
        }
        case fixedHeatTransferCoeff:
        {
            writeEntry(os, "h", h_);
            writeEntry(os, Ta_());

            if (relaxation_ < 1)
            {
                writeEntry(os, "relaxation", relaxation_);
            }

            if (emissivity_ > 0)
            {
                writeEntry(os, "emissivity", emissivity_);
            }

            if (thicknessLayers_.size())
            {
                writeEntry(os, "thicknessLayers", thicknessLayers_);
                writeEntry(os, "kappaLayers", kappaLayers_);
            }

            break;
        }
    }

    writeEntry(os, "qr", qrName_);

    if (qrName_ != "none")
    {
        writeEntry(os, "qrRelaxation", qrRelaxation_);

        writeEntry(os, "qrPrevious", qrPrevious_);
    }

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "refGradient", refGrad());
    writeEntry(os, "valueFraction", valueFraction());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousMediaWallHeatTransferFvPatchScalarField
    );
}

// ************************************************************************* //
