/*---------------------------------------------------------------------------*\
    Solver-side porous-media closure input validation.
    Dictionary keywords are printed together with manuscript symbols and units.
\*---------------------------------------------------------------------------*/
#ifndef porousClosureInputValidation_H
#define porousClosureInputValidation_H

#include "wordList.H"
#include "dictionary.H"
#include "IStringStream.H"
#include "OStringStream.H"
#include <sstream>
#include <iomanip>
#include <vector>
#include "DynamicList.H"
#include "scalar.H"
#include "error.H"
#include <cmath>

namespace Foam
{

inline void requirePorousModelInputs
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

inline word porousModelSelection
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

inline void checkPorousScalarInput
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

inline void validatePorousCommonInputs(const dictionary& dict)
{
    const string context("common porous-medium properties");
    requirePorousModelInputs
    (
        dict, context,
        "(Porosity ParticleDiameter SolidDensity SolidHeatCapacity)",
        "Porosity = epsilon [-]; ParticleDiameter = d_p [m]; "
        "SolidDensity = rho_s [kg/m3]; SolidHeatCapacity = c_p,s [J/kg/K]"
    );
    checkPorousScalarInput(dict, "Porosity", context, "epsilon [-]");
    if (readScalar(dict.lookup("Porosity")) >= 1)
    {
        FatalIOErrorInFunction(dict)
            << "Porosity (epsilon) must satisfy 0 < epsilon < 1."
            << exit(FatalIOError);
    }
    checkPorousScalarInput(dict, "ParticleDiameter", context, "d_p [m]");
    checkPorousScalarInput(dict, "SolidDensity", context, "rho_s [kg/m3]");
    checkPorousScalarInput(dict, "SolidHeatCapacity", context, "c_p,s [J/kg/K]");
}

inline void validatePorousMomentumInputs(const dictionary& dict)
{
    const word model = porousModelSelection
    (dict, "pressureDropClosureModel", "(None Darcy Forchheimer Ergun)");
    const string context("pressureDropClosureModel = " + model);
    if (model == "Darcy" || model == "Forchheimer")
    {
        requirePorousModelInputs
        (
            dict, context,
            model == "Darcy" ? "(DarcyCoeff)" : "(DarcyCoeff ForchCoeff)",
            "DarcyCoeff = D [1/m2]; ForchCoeff = F [1/m] (Forchheimer only)"
        );
        checkPorousScalarInput(dict, "DarcyCoeff", context, "D [1/m2]", true);
        if (model == "Forchheimer")
            checkPorousScalarInput(dict, "ForchCoeff", context, "F [1/m]", true);
    }
}

inline void validatePorousThermalInputs(const dictionary& dict)
{
    const word model = porousModelSelection
    (
        dict, "ThermalConductivityClosureModel",
        "(Const ConstAnisotropic Standard Bunnell Yagi Wakao)"
    );
    const string context("ThermalConductivityClosureModel = " + model);
    if (model == "Const")
    {
        requirePorousModelInputs(dict, context, "(EffectiveThermalConductivity)",
            "EffectiveThermalConductivity = k_eff [W/m/K]");
        checkPorousScalarInput(dict, "EffectiveThermalConductivity", context,
            "k_eff [W/m/K]");
    }
    else if (model == "ConstAnisotropic")
    {
        requirePorousModelInputs(dict, context,
            "(EffectiveRadialThermalConductivity EffectiveAxialThermalConductivity AxialDirection)",
            "EffectiveRadialThermalConductivity = k_eff,r [W/m/K]; "
            "EffectiveAxialThermalConductivity = k_eff,z [W/m/K]; "
            "AxialDirection = e_z, finite nonzero vector [-]");
        checkPorousScalarInput(dict, "EffectiveRadialThermalConductivity", context,
            "k_eff,r [W/m/K]");
        checkPorousScalarInput(dict, "EffectiveAxialThermalConductivity", context,
            "k_eff,z [W/m/K]");
        // Vector parsing, validation and normalization are performed by the reader.
    }
    else
    {
        const char* names = "(SolidThermalConductivity)";
        const char* description = "SolidThermalConductivity = k_s [W/m/K]";
        if (model == "Bunnell")
        {
            names = "(SuperficialVelocity)";
            description = "SuperficialVelocity = |<u>| [m/s]";
        }
        else if (model == "Yagi")
        {
            names = "(SolidThermalConductivity SuperficialVelocity YagiAlphaBeta YagiPhi)";
            description = "SolidThermalConductivity = k_s [W/m/K]; "
                "SuperficialVelocity = |<u>| [m/s]; "
                "YagiAlphaBeta = (alpha beta)_Y [-]; YagiPhi = phi_Y [-]";
        }
        else if (model == "Wakao")
        {
            names = "(SolidThermalConductivity SuperficialVelocity StagnantThermalConductivityClosureModel)";
            description = "SolidThermalConductivity = k_s [W/m/K]; "
                "SuperficialVelocity = |<u>| [m/s]; "
                "StagnantThermalConductivityClosureModel = Kunii, Krupiczka or ZS "
                "for k_eff^0; optional WakaoMixCoeff = C_mix [-], default 0.1";
        }
        requirePorousModelInputs(dict, context, names, description);
        if (model != "Bunnell")
            checkPorousScalarInput(dict, "SolidThermalConductivity", context, "k_s [W/m/K]");
        if (model != "Standard")
            checkPorousScalarInput(dict, "SuperficialVelocity", context, "|<u>| [m/s]", true);
        if (model == "Yagi")
        {
            checkPorousScalarInput(dict, "YagiAlphaBeta", context, "(alpha beta)_Y [-]", true);
            checkPorousScalarInput(dict, "YagiPhi", context, "phi_Y [-]");
        }
        if (model == "Wakao")
        {
            const word stagnant = porousModelSelection
            (dict, "StagnantThermalConductivityClosureModel", "(Kunii Krupiczka ZS)");
            if (stagnant == "Kunii")
            {
                const string kunii("ThermalConductivityClosureModel = Wakao, stagnant model = Kunii");
                requirePorousModelInputs(dict, kunii, "(KuniiPhi KuniiGamma)",
                    "KuniiPhi = phi_K [-]; KuniiGamma = gamma_K [-]");
                checkPorousScalarInput(dict, "KuniiPhi", kunii, "phi_K [-]");
                checkPorousScalarInput(dict, "KuniiGamma", kunii, "gamma_K [-]", true);
            }
            if (dict.found("WakaoMixCoeff"))
                checkPorousScalarInput(dict, "WakaoMixCoeff", context, "C_mix [-]", true);
        }
    }
}

} // End namespace Foam
#endif
