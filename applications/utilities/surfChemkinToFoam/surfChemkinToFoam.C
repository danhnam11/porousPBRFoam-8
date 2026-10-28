/* surfChemkinToFoam-new: checked, backward-compatible constant-thermo
   surface mechanism conversion for OpenFOAM Foundation 8. */
#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "atomicWeights.H"
#include "surfaceMechanism.H"
#include <cstdio>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>

using namespace Foam;
namespace sc = surfaceConversion;

namespace
{
std::string canonical(const fileName& file)
{
    char* existing=::realpath(file.c_str(),0);
    if(existing) {std::string result(existing);std::free(existing);return result;}
    const fileName parent=file.path().empty()?fileName("."):file.path();
    char* directory=::realpath(parent.c_str(),0);
    if(!directory)throw std::runtime_error("Output directory does not exist: "+std::string(parent.c_str()));
    std::string result(directory);std::free(directory);return result+"/"+file.name().c_str();
}
double property(const dictionary& dict,const char* key,double fallback)
{
    return dict.lookupOrDefault<scalar>(key,fallback);
}
std::string chemicalFormula(std::string name)
{
    if(name.size()>=3 && (name.substr(name.size()-3)=="(s)" || name.substr(name.size()-3)=="(S)"))name.erase(name.size()-3);
    // Preserve chemical formulas such as CO; only normalize known bare metals.
    if(name=="NI")return "Ni";
    if(name=="PT")return "Pt";
    if(name=="RH")return "Rh";
    if(name=="RU")return "Ru";
    return name;
}
double molecularWeight(const std::string& name)
{
    const std::string formula=chemicalFormula(name);double weight=0;
    for(std::size_t i=0;i<formula.size();)
    {
        if(!std::isupper(static_cast<unsigned char>(formula[i])))throw std::runtime_error("Cannot infer molecular weight for "+name+"; specify molWeight in surfaceThermoOverrides.");
        std::string element(1,formula[i++]);
        if(i<formula.size() && std::islower(static_cast<unsigned char>(formula[i])))element+=formula[i++];
        std::string n;while(i<formula.size() && std::isdigit(static_cast<unsigned char>(formula[i])))n+=formula[i++];
        const double count=n.empty()?1:sc::number(n,"species name",0);
        const word e(element.c_str());
        if(!atomicWeights.found(e))throw std::runtime_error("Unknown element in species name "+name+"; specify molWeight in surfaceThermoOverrides.");
        weight+=count*atomicWeights[e];
    }
    if(weight<=0)throw std::runtime_error("Nonpositive molecular weight: "+name);
    return weight;
}
void validProperties(const std::string& name,double rho,double mu,double Pr,double Cp,double Cv,double Hf,double mw)
{
    if(!std::isfinite(rho)||!std::isfinite(mu)||!std::isfinite(Pr)||!std::isfinite(Cp)||!std::isfinite(Cv)||!std::isfinite(Hf)||!std::isfinite(mw)
      ||rho<=0||mu<0||Pr<=0||Cp<0||Cv<0||mw<=0)
        throw std::runtime_error("Invalid constant thermophysical property for "+name+" (rho, Pr, molWeight > 0; mu, Cp, Cv >= 0).");
}
std::string thermoOutput(const sc::Mechanism& m,const std::vector<sc::ThermoRecord>& records,const dictionary& config,const std::string& object)
{
    const dictionary& defaults=config.subDict("surfaceThermoDefaults");
    // Retain the original mandatory rho/mu/Pr entries and optional Cp/Cv/Hf.
    const double rho0=readScalar(defaults.lookup("rho"));
    const double mu0=readScalar(defaults.lookup("mu"));
    const double Pr0=readScalar(defaults.lookup("Pr"));
    const double cp0=property(defaults,"Cp",0),cv0=property(defaults,"Cv",0),hf0=property(defaults,"Hf",0);
    const dictionary* overrides=config.found("surfaceThermoOverrides")?&config.subDict("surfaceThermoOverrides"):0;
    if(overrides)
    {
        const wordList names=overrides->toc();
        forAll(names,i)
        {
            if(std::find(m.species.begin(),m.species.end(),std::string(names[i].c_str()))==m.species.end())
                throw std::runtime_error("surfaceThermoOverrides species is not in the mapped SITE list: "+std::string(names[i].c_str()));
            const dictionary& values=overrides->subDict(names[i]);const wordList keys=values.toc();
            forAll(keys,j)
            {
                if(keys[j]!="rho" && keys[j]!="mu" && keys[j]!="Pr" && keys[j]!="Cp" && keys[j]!="Cv" && keys[j]!="Hf" && keys[j]!="molWeight")
                    throw std::runtime_error("Unknown surfaceThermoOverrides property: "+std::string(keys[j].c_str()));
            }
        }
    }
    std::set<std::string> available;for(std::size_t i=0;i<records.size();++i)available.insert(records[i].name);
    std::ostringstream out;out<<std::setprecision(10)<<sc::header(object)<<"species\n"<<m.species.size()<<"\n(\n";
    for(std::size_t i=0;i<m.species.size();++i)out<<"    "<<m.species[i]<<"\n";
    out<<");\n\n";
    for(std::size_t i=0;i<m.species.size();++i)
    {
        const std::string& name=m.species[i];
        if(!available.count(name))throw std::runtime_error("Missing thermodynamic record for SITE species: "+name);
        double rho=rho0,mu=mu0,Pr=Pr0,Cp=cp0,Cv=cv0,Hf=hf0,mw=0;
        bool givenWeight=false;
        if(overrides && overrides->found(name.c_str()))
        {
            const dictionary& values=overrides->subDict(name.c_str());
            rho=property(values,"rho",rho);mu=property(values,"mu",mu);Pr=property(values,"Pr",Pr);
            Cp=property(values,"Cp",Cp);Cv=property(values,"Cv",Cv);Hf=property(values,"Hf",Hf);
            if(values.found("molWeight")){mw=readScalar(values.lookup("molWeight"));givenWeight=true;}
        }
        if(!givenWeight)mw=molecularWeight(name);
        validProperties(name,rho,mu,Pr,Cp,Cv,Hf,mw);
        out<<name<<"\n{\n    specie\n    {\n        molWeight "<<mw<<";\n        siteDensity "<<m.siteDensity<<";\n    }\n"
           <<"    equationOfState\n    {\n        rho "<<rho<<";\n    }\n"
           <<"    thermodynamics\n    {\n        Cp "<<Cp<<";\n        Hf "<<Hf<<";\n        Cv "<<Cv<<";\n    }\n"
           <<"    transport\n    {\n        mu "<<mu<<";\n        Pr "<<Pr<<";\n    }\n}\n\n";
    }
    return out.str();
}
std::string stage(const std::string& destination,const std::string& contents)
{
    std::string pattern=destination+".tmp.XXXXXX";std::vector<char> name(pattern.begin(),pattern.end());name.push_back(0);
    int fd=::mkstemp(&name[0]);if(fd<0)throw std::runtime_error("Cannot create output beside "+destination);
    const std::string temp(&name[0]);std::size_t offset=0;
    while(offset<contents.size())
    {
        const ssize_t count=::write(fd,contents.data()+offset,contents.size()-offset);
        if(count<0 && errno==EINTR)continue;
        if(count<=0){::close(fd);std::remove(temp.c_str());throw std::runtime_error("Write error: "+destination);}
        offset+=count;
    }
    if(::close(fd)!=0){std::remove(temp.c_str());throw std::runtime_error("Close error: "+destination);}
    return temp;
}
void writePair(const std::string& a,const std::string& at,const std::string& b,const std::string& bt)
{
    const std::string suffix=".backup."+std::to_string(static_cast<long long>(::getpid()));
    std::string ta,tb;bool savedA=false,savedB=false,wroteA=false,wroteB=false;
    try
    {
        ta=stage(a,at);tb=stage(b,bt);struct stat st;
        if(::stat(a.c_str(),&st)==0)
        {
            if(!S_ISREG(st.st_mode))throw std::runtime_error("Output is not a regular file: "+a);
            if(::access((a+suffix).c_str(),F_OK)==0)throw std::runtime_error("Backup path already exists: "+a+suffix);
            if(std::rename(a.c_str(),(a+suffix).c_str()))throw std::runtime_error("Cannot preserve previous output: "+a);
            savedA=true;
        }
        if(::stat(b.c_str(),&st)==0)
        {
            if(!S_ISREG(st.st_mode))throw std::runtime_error("Output is not a regular file: "+b);
            if(::access((b+suffix).c_str(),F_OK)==0)throw std::runtime_error("Backup path already exists: "+b+suffix);
            if(std::rename(b.c_str(),(b+suffix).c_str()))throw std::runtime_error("Cannot preserve previous output: "+b);
            savedB=true;
        }
        if(std::rename(ta.c_str(),a.c_str()))throw std::runtime_error("Cannot install output: "+a);
        wroteA=true;
        if(std::rename(tb.c_str(),b.c_str()))throw std::runtime_error("Cannot install output: "+b);
        wroteB=true;
        if(savedA)std::remove((a+suffix).c_str());
        if(savedB)std::remove((b+suffix).c_str());
    }
    catch(...)
    {
        if(wroteA)std::remove(a.c_str());
        if(wroteB)std::remove(b.c_str());
        if(savedA)std::rename((a+suffix).c_str(),a.c_str());
        if(savedB)std::rename((b+suffix).c_str(),b.c_str());
        if(!ta.empty())std::remove(ta.c_str());
        if(!tb.empty())std::remove(tb.c_str());
        throw;
    }
}
}
int main(int argc,char* argv[])
{
    #include "removeCaseOptions.H"
    argList::validArgs.append("CHEMKIN surface chemistry file");
    argList::validArgs.append("CHEMKIN surface thermodynamics file");
    argList::validArgs.append("surfaceThermoDict");
    argList::validArgs.append("OpenFOAM surface chemistry file");
    argList::validArgs.append("OpenFOAM surface thermodynamics file");
    argList::addBoolOption("checkOnly","validate both inputs without writing outputs");
    argList::addOption("siteSpecies","name","map bare * to this name and X(*) to X(s)");
    argList::addBoolOption("newFormat","legacy compatibility option; NASA7 formats are detected automatically");
    argList args(argc,argv);
    try
    {
        sc::Settings settings;
        if(args.optionFound("siteSpecies"))
        { settings.bareSite=args.optionLookupOrDefault<word>("siteSpecies",word::null).c_str();if(!sc::safeName(settings.bareSite))throw std::runtime_error("Invalid -siteSpecies name"); }
        IFstream input(args[3]);if(!input.good())throw std::runtime_error("Cannot open surfaceThermoDict: "+std::string(args[3].c_str()));
        dictionary config(input);
        if(config.found("speciesMap"))
        {
            const dictionary& map=config.subDict("speciesMap");const wordList keys=map.toc();
            forAll(keys,i)
            {
                Foam::string value;map.lookup(keys[i])>>value;
                if(!sc::safeName(value.c_str()))throw std::runtime_error("Invalid speciesMap target: "+std::string(value.c_str()));
                settings.aliases[keys[i].c_str()]=value.c_str();
            }
        }
        sc::Mechanism mechanism=sc::readMechanism(args[1].c_str(),settings);
        const std::vector<sc::ThermoRecord> thermo=sc::readThermo(args[2].c_str(),settings,mechanism);
        const std::string a=canonical(fileName(args[4])),b=canonical(fileName(args[5]));
        if(a==b)throw std::runtime_error("Reaction and thermo outputs must be different files");
        for(int i=1;i<=3;++i)if(a==canonical(fileName(args[i])) || b==canonical(fileName(args[i])))throw std::runtime_error("An output path aliases an input file");
        const std::string reactionText=sc::reactionOutput(mechanism,fileName(args[4]).name().c_str());
        const std::string thermoText=thermoOutput(mechanism,thermo,config,fileName(args[5]).name().c_str());
        std::size_t stick=0,coverage=0,pkc=0;
        for(std::size_t i=0;i<mechanism.reactions.size();++i)
        {stick+=mechanism.reactions[i].stick;coverage+=!mechanism.reactions[i].coverage.empty();pkc+=mechanism.reactions[i].pkc;}
        Info<<"Surface mechanism validation passed."<<nl
            <<"  SITE species: "<<mechanism.species.size()<<nl
            <<"  Reactions: "<<mechanism.reactions.size()<<" (STICK "<<stick<<", COV "<<coverage<<", PKC "<<pkc<<")"<<nl
            <<"  Energy unit: "<<mechanism.energyUnit.c_str()<<nl
            <<"  Amount unit: "<<(mechanism.molecules?"MOLECULES":"MOLES")<<nl
            <<"  Site density [kmol/m2]: "<<mechanism.siteDensity<<nl
            <<"  Thermo output: constant properties from surfaceThermoDefaults and surfaceThermoOverrides."<<nl
            <<"  NASA7 records are checked; NASA coefficients are not used to replace Cp or Hf."<<nl;
        if(args.optionFound("newFormat"))Info<<"Note: -newFormat is retained for command-line compatibility; NASA7 records are detected automatically."<<nl;
        for(std::size_t i=0;i<mechanism.warnings.size();++i)Info<<"Warning: "<<mechanism.warnings[i].c_str()<<nl;
        for(std::map<std::string,std::string>::const_iterator it=settings.aliases.begin();it!=settings.aliases.end();++it)Info<<"  Species mapping: "<<it->first.c_str()<<" -> "<<it->second.c_str()<<nl;
        if(!settings.bareSite.empty())Info<<"  Bare-site mapping: * -> "<<settings.bareSite.c_str()<<"; X(*) -> X(s)"<<nl;
        if(args.optionFound("checkOnly"))Info<<"Check-only mode: no output files written."<<nl;
        else {writePair(a,reactionText,b,thermoText);Info<<"Written: "<<args[4]<<" and "<<args[5]<<nl;}
        return 0;
    }
    catch(const std::exception& error)
    {Info<<"Conversion failed: "<<error.what()<<nl;return 1;}
}
