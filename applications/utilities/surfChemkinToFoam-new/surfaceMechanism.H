// Shared CHEMKIN surface parser. C++11; no lexer or external runtime required.
#ifndef surfaceMechanism_H
#define surfaceMechanism_H
#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace surfaceConversion
{
struct Coverage { std::string species; double eta, mu, epsilon; };
struct Reaction
{
    std::string equation;
    double A, beta, Ta, gasOrder, surfaceOrder;
    bool stick, pkc;
    std::vector<Coverage> coverage;
    std::size_t line;
};
struct Mechanism
{
    std::vector<std::string> species;
    std::vector<Reaction> reactions;
    double siteDensity;
    double energyToKelvin;
    bool molecules;
    std::string energyUnit;
    std::vector<std::string> warnings;
};
struct ThermoRecord
{
    std::string name;
    std::vector<double> coefficients;
};
struct Settings
{
    std::map<std::string,std::string> aliases;
    std::string bareSite;
};
inline std::string trim(const std::string& s)
{
    const std::size_t a=s.find_first_not_of(" \t\r\n");
    return a==std::string::npos ? "" : s.substr(a,s.find_last_not_of(" \t\r\n")-a+1);
}
inline std::string upper(std::string s)
{
    for(std::size_t i=0;i<s.size();++i) s[i]=std::toupper(static_cast<unsigned char>(s[i]));
    return s;
}
inline std::string clean(const std::string& s) { return trim(s.substr(0,s.find('!'))); }
inline std::vector<std::string> tokens(const std::string& s)
{
    std::istringstream in(s); std::vector<std::string> result; std::string t;
    while(in>>t) result.push_back(t);
    return result;
}
inline void fail(const std::string& file,std::size_t line,const std::string& why)
{
    std::ostringstream out; out<<file<<":"<<line<<": "<<why;
    throw std::runtime_error(out.str());
}
inline double number(std::string s,const std::string& file,std::size_t line)
{
    std::replace(s.begin(),s.end(),'D','E'); std::replace(s.begin(),s.end(),'d','e');
    char* end=0; errno=0; const double value=std::strtod(s.c_str(),&end);
    if(s.empty() || end==s.c_str() || *end || errno==ERANGE || !std::isfinite(value))
        fail(file,line,"invalid number '"+s+"'");
    return value;
}
inline std::string mapped(const std::string& s,const Settings& settings)
{
    std::map<std::string,std::string>::const_iterator it=settings.aliases.find(s);
    if(it!=settings.aliases.end()) return it->second;
    if(!settings.bareSite.empty())
    {
        if(s=="*") return settings.bareSite;
        if(s.size()>3 && s.substr(s.size()-3)=="(*)") return s.substr(0,s.size()-3)+"(s)";
    }
    return s;
}
inline bool safeName(const std::string& s)
{
    if(s.empty()) return false;
    const std::string invalid=" \t\r\n\"{};/+*=<>";
    return s.find_first_of(invalid)==std::string::npos;
}
inline std::string readLines(const std::string& file,std::vector<std::string>& lines)
{
    std::ifstream in(file.c_str());
    if(!in) throw std::runtime_error("Cannot open input file: "+file);
    std::string l; while(std::getline(in,l)) lines.push_back(l);
    if(in.bad()) throw std::runtime_error("Read error: "+file);
    return file;
}
inline std::string siteToken(std::string token,const std::string& file,std::size_t line)
{
    if(!token.empty() && (token.back()==',' || token.back()==';')) token.pop_back();
    const std::size_t slash=token.find('/');
    if(slash!=std::string::npos)
    {
        if(token.back()!='/' || slash==0) fail(file,line,"invalid SITE species/occupancy: "+token);
        if(number(token.substr(slash+1,token.size()-slash-2),file,line)!=1)
            fail(file,line,"multi-site occupancy is not supported by the current output model");
        token=token.substr(0,slash);
    }
    return token;
}
struct Side { std::string text; double gas, surface; };
inline Side side(const std::string& text,const Mechanism& m,const Settings& settings,const std::string& file,std::size_t line)
{
    Side out; out.gas=out.surface=0;
    std::size_t start=0;
    while(start<text.size())
    {
        const std::size_t end=text.find('+',start);
        std::string term=trim(text.substr(start,end==std::string::npos?end:end-start));
        if(term.empty()) fail(file,line,"empty reaction term");
        double coefficient=1;
        if(std::isdigit(static_cast<unsigned char>(term[0])) || term[0]=='.')
        {
            std::size_t i=0;
            while(i<term.size() && (std::isdigit(static_cast<unsigned char>(term[i])) || term[i]=='.')) ++i;
            coefficient=number(term.substr(0,i),file,line); term=trim(term.substr(i));
        }
        if(coefficient<=0 || term.empty() || tokens(term).size()!=1) fail(file,line,"invalid stoichiometric term");
        const std::string name=mapped(term,settings);
        const bool isSurface=std::find(m.species.begin(),m.species.end(),name)!=m.species.end();
        if(!safeName(name)) fail(file,line,"invalid/unmapped species name: "+name);
        if(!isSurface && (name.find("(s)")!=std::string::npos || name.find("(S)")!=std::string::npos || name.find("(*)")!=std::string::npos))
            fail(file,line,"surface species not declared in SITE: "+name);
        if(isSurface) out.surface+=coefficient; else out.gas+=coefficient;
        if(!out.text.empty()) out.text+=" + ";
        if(coefficient!=1) { std::ostringstream v;v<<std::setprecision(10)<<coefficient;out.text+=v.str(); }
        out.text+=name;
        if(end==std::string::npos) break;
        start=end+1;
        if(start==text.size()) fail(file,line,"trailing '+' in reaction");
    }
    if(out.text.empty()) fail(file,line,"empty reaction side");
    return out;
}
inline Mechanism readMechanism(const std::string& file,const Settings& settings)
{
    Mechanism m; m.siteDensity=0; m.molecules=false; m.energyToKelvin=1000.0/8.31446261815324;
    m.energyUnit="KJOULES/MOLE";
    std::vector<std::string> lines;readLines(file,lines);
    bool inSite=false,siteEnded=false,inReactions=false,reactionsEnded=false,seenSite=false;
    std::size_t pending=static_cast<std::size_t>(-1);
    for(std::size_t i=0;i<lines.size();++i)
    {
        const std::string text=clean(lines[i]);if(text.empty())continue;
        const std::string u=upper(text);
        if(u.compare(0,4,"SITE")==0)
        {
            if(seenSite)fail(file,i+1,"multiple SITE phases are not supported by the current solver");
            seenSite=true;inSite=true;
            const std::size_t pos=u.find("SDEN");
            if(pos==std::string::npos) fail(file,i+1,"SITE header must contain SDEN/value/");
            const std::size_t a=text.find('/',pos),b=a==std::string::npos?a:text.find('/',a+1);
            if(a==std::string::npos || b==std::string::npos)fail(file,i+1,"malformed SDEN");
            m.siteDensity=10*number(trim(text.substr(a+1,b-a-1)),file,i+1);
            if(m.siteDensity<=0)fail(file,i+1,"SDEN must be positive");
            if(!trim(text.substr(b+1)).empty())fail(file,i+1,"put SITE species on the lines following the SITE header");
            continue;
        }
        if(inSite)
        {
            const std::vector<std::string> names=tokens(text);
            for(std::size_t k=0;k<names.size();++k)
            {
                if(upper(names[k])=="END") {inSite=false;siteEnded=true;break;}
                if(names[k].size()>2 && names[k][0]=='/' && names[k].back()=='/')
                {
                    if(number(names[k].substr(1,names[k].size()-2),file,i+1)!=1)fail(file,i+1,"non-unit site occupancy is unsupported");
                    continue;
                }
                const std::string name=mapped(siteToken(names[k],file,i+1),settings);
                if(!safeName(name))fail(file,i+1,"invalid SITE name; supply speciesMap or -siteSpecies: "+name);
                if(std::find(m.species.begin(),m.species.end(),name)!=m.species.end())fail(file,i+1,"duplicate/mapping-colliding SITE species: "+name);
                m.species.push_back(name);
            }
            continue;
        }
        if(u.compare(0,9,"REACTIONS")==0)
        {
            if(inReactions || reactionsEnded) fail(file,i+1,"multiple REACTIONS blocks are unsupported");
            if(!siteEnded || m.species.empty())fail(file,i+1,"a complete SITE block is required before REACTIONS");
            inReactions=true;bool specifiedEnergy=false,specifiedAmount=false;
            const std::vector<std::string> units=tokens(u);
            for(std::size_t k=1;k<units.size();++k)
            {
                const std::string& unit=units[k]; double factor=0;
                if(unit=="CAL/MOLE" || unit=="CAL/MOL")factor=4.184/8.31446261815324;
                else if(unit=="KCAL/MOLE" || unit=="KCAL/MOL")factor=4184/8.31446261815324;
                else if(unit=="JOULES/MOLE" || unit=="JOULES/MOL")factor=1/8.31446261815324;
                else if(unit=="KJOULES/MOLE" || unit=="KJOULES/MOL")factor=1000/8.31446261815324;
                else if(unit=="KELVINS")factor=1;
                else if(unit=="EVOLTS")factor=96485.33212331002/8.31446261815324;
                if(factor)
                {
                    if(specifiedEnergy)fail(file,i+1,"multiple activation-energy units");
                    specifiedEnergy=true;m.energyUnit=unit;m.energyToKelvin=factor;
                }
                else if(unit=="MOLES" || unit=="MOLECULES")
                {
                    if(specifiedAmount)fail(file,i+1,"multiple amount units");
                    specifiedAmount=true;m.molecules=unit=="MOLECULES";
                }
                else if(unit!="MWOFF")fail(file,i+1,"unsupported REACTIONS option: "+unit);
            }
            if(!specifiedEnergy)m.warnings.push_back("No energy unit specified: retaining legacy KJOULES/MOLE assumption; specify it explicitly.");
            continue;
        }
        if(!inReactions)
        {
            if(u.compare(0,4,"BULK")==0)fail(file,i+1,"BULK phases are unsupported");
            continue; // embedded thermo data are read from the explicit thermo input
        }
        if(u=="END") {reactionsEnded=true;inReactions=false;continue;}
        if(reactionsEnded)fail(file,i+1,"unexpected content after REACTIONS END");
        if(u=="STICK" || u=="PKC" || u.compare(0,3,"COV")==0)
        {
            if(pending==static_cast<std::size_t>(-1))fail(file,i+1,"modifier without a preceding reaction");
            Reaction& r=m.reactions[pending];
            if(u=="STICK") {if(r.stick||r.pkc)fail(file,i+1,"duplicate or conflicting rate modifier");r.stick=true;}
            else if(u=="PKC") {if(r.pkc||r.stick)fail(file,i+1,"duplicate or conflicting rate modifier");r.pkc=true;}
            else
            {
                std::string v=text;std::replace(v.begin(),v.end(),'/',' ');
                std::vector<std::string> t=tokens(v);
                if(t.size()!=5)fail(file,i+1,"COV requires species, eta, mu and epsilon");
                Coverage c;c.species=mapped(t[1],settings);
                if(std::find(m.species.begin(),m.species.end(),c.species)==m.species.end())fail(file,i+1,"COV species is not in SITE: "+c.species);
                for(std::size_t j=0;j<r.coverage.size();++j)if(r.coverage[j].species==c.species)fail(file,i+1,"duplicate COV species");
                c.eta=number(t[2],file,i+1);c.mu=number(t[3],file,i+1);c.epsilon=number(t[4],file,i+1)*m.energyToKelvin;
                r.coverage.push_back(c);
            }
            continue;
        }
        const std::vector<std::string> t=tokens(text);
        if(t.size()<4)fail(file,i+1,"unsupported modifier or malformed reaction: "+text);
        Reaction r;r.A=number(t[t.size()-3],file,i+1);r.beta=number(t[t.size()-2],file,i+1);
        r.Ta=number(t.back(),file,i+1)*m.energyToKelvin;r.stick=r.pkc=false;r.line=i+1;
        if(r.A<0)fail(file,i+1,"negative Arrhenius prefactor");
        std::string equation;for(std::size_t k=0;k+3<t.size();++k){if(k)equation+=' ';equation+=t[k];}
        const std::size_t arrow=equation.find("=>");
        if(arrow==std::string::npos || equation.find("<=>")!=std::string::npos)
            fail(file,i+1,"only explicit irreversible '=>' reactions are supported; reversible reactions must not be silently converted");
        const Side lhs=side(equation.substr(0,arrow),m,settings,file,i+1);
        const Side rhs=side(equation.substr(arrow+2),m,settings,file,i+1);
        if(std::fabs(lhs.surface-rhs.surface)>1e-9)fail(file,i+1,"reaction does not conserve the number of surface sites");
        r.gasOrder=lhs.gas;r.surfaceOrder=lhs.surface;r.equation=lhs.text+" = "+rhs.text;
        m.reactions.push_back(r);pending=m.reactions.size()-1;
    }
    if(!seenSite || !siteEnded || !reactionsEnded || m.reactions.empty())fail(file,lines.size(),"missing/empty SITE or REACTIONS block, or missing END");
    for(std::size_t i=0;i<m.reactions.size();++i)
    {
        Reaction& r=m.reactions[i];
        if(r.pkc && !r.coverage.empty())fail(file,r.line,"PKC with COV is unsupported");
        if(r.stick && !r.coverage.empty())fail(file,r.line,"combined STICK and COV requires a combined library rate model; unsupported");
        if(!r.stick)
        {
            r.A*=10*std::pow(1e-3,r.gasOrder)*std::pow(1e-1,r.surfaceOrder);
            if(m.molecules)r.A*=std::pow(6.02214076e23,r.gasOrder+r.surfaceOrder-1);
        }
        if(!std::isfinite(r.A))fail(file,r.line,"prefactor unit conversion overflow");
    }
    return m;
}
inline std::vector<ThermoRecord> readThermo(const std::string& file,const Settings& settings,Mechanism& mechanism)
{
    std::vector<std::string> lines;readLines(file,lines);std::vector<ThermoRecord> out;
    bool active=false,ended=false;std::set<std::string> seen;
    for(std::size_t i=0;i<lines.size();++i)
    {
        const std::string text=clean(lines[i]);if(text.empty())continue;
        const std::string u=upper(text);
        if(!active) {if(u=="THERMO" || u=="THERMO ALL")active=true;continue;}
        if(u=="END"){ended=true;break;}
        const std::vector<std::string> t=tokens(text);
        if(t.size()==3 && !t[0].empty() && (std::isdigit(static_cast<unsigned char>(t[0][0]))||t[0][0]=='.'))
        {for(std::size_t k=0;k<3;++k)number(t[k],file,i+1);continue;}
        ThermoRecord e;
        const std::vector<std::string> header=tokens(lines[i].substr(0,18));
        if(header.empty())fail(file,i+1,"missing NASA7 species name");
        e.name=mapped(header[0],settings);
        const bool required=std::find(mechanism.species.begin(),mechanism.species.end(),e.name)!=mechanism.species.end();
        if(!seen.insert(e.name).second)
        {
            if(required)fail(file,i+1,"duplicate/mapping-colliding thermo species: "+e.name);
            mechanism.warnings.push_back("Ignoring repeated thermo record not present in SITE: "+e.name);
        }
        for(int row=0;row<3;++row)
        {
            if(++i>=lines.size())fail(file,i,"truncated NASA7 coefficient record: "+e.name);
            std::string coefficient=lines[i].substr(0,75);
            std::replace(coefficient.begin(),coefficient.end(),'D','E');std::replace(coefficient.begin(),coefficient.end(),'d','e');
            const char* cursor=coefficient.c_str();int count=0;
            while(*cursor)
            {
                while(*cursor && std::isspace(static_cast<unsigned char>(*cursor)))++cursor;
                if(!*cursor)break;
                char* end=0;errno=0;double v=std::strtod(cursor,&end);
                if(end==cursor || errno==ERANGE || !std::isfinite(v))fail(file,i+1,"invalid NASA7 coefficient for "+e.name);
                e.coefficients.push_back(v);++count;cursor=end;
            }
            if(count!=(row==2?4:5))fail(file,i+1,"expected NASA7 coefficient row with "+std::string(row==2?"4":"5")+" values for "+e.name);
        }
        out.push_back(e);
    }
    if(!active || out.empty())fail(file,lines.size(),"missing/empty THERMO block");
    if(!ended)mechanism.warnings.push_back("THERMO block ended at EOF after a complete NASA7 record (legacy compatibility); add END to the input.");
    return out;
}
inline std::string header(const std::string& object)
{
    return "/* Generated by surfChemkinToFoam-new */\nFoamFile\n{\n    version 2.0;\n    format ascii;\n    class dictionary;\n    location \"constant\";\n    object \""+object+"\";\n}\n\n";
}
inline std::string reactionOutput(const Mechanism& m,const std::string& object)
{
    std::ostringstream out;out<<std::setprecision(10)<<header(object)<<"reactions\n{\n";
    for(std::size_t i=0;i<m.reactions.size();++i)
    {
        const Reaction& r=m.reactions[i];
        const std::string type=!r.coverage.empty()?"CoverageModificationArrhenius":r.pkc?"Takahashi":r.stick?"StickArrhenius":"Arrhenius";
        out<<"    un-named-reaction-"<<i<<"\n    {\n        type irreversibleSurface"<<type<<"Surface;\n"
           <<"        reaction \""<<r.equation<<"\";\n        A "<<r.A<<";\n        beta "<<r.beta<<";\n        Ta "<<r.Ta<<";\n";
        if(!r.coverage.empty())
        {
            out<<"        coeffs\n        "<<m.species.size()<<"\n        (\n";
            for(std::size_t k=0;k<m.species.size();++k)
            {
                Coverage c;c.eta=c.mu=c.epsilon=0;c.species=m.species[k];
                for(std::size_t j=0;j<r.coverage.size();++j)if(r.coverage[j].species==c.species)c=r.coverage[j];
                out<<"            ("<<c.species<<" ("<<c.eta<<" "<<c.mu<<" "<<c.epsilon<<"))\n";
            }
            out<<"        );\n";
        }
        out<<"    }\n";
    }
    out<<"}\n\nTlow 0;\nThigh 1e300;\n";return out.str();
}
}
#endif
