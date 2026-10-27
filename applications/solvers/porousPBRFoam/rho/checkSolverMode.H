{
            const IOdictionary porousProperties
            (
                IOobject
                (
                    "porousMediaProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            const word mode =
                porousProperties.lookupOrDefault<word>("solverMode", "PM");
            if (mode != "PM")
            {
                FatalErrorInFunction
                    << "heRhoThermo supports solverMode PM only. "
                    << "Use hePsiThermo for PMSCHEM surface chemistry."
                    << exit(FatalError);
            }
        }
