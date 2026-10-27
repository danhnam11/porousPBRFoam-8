word thermoType;
    {
        const IOdictionary properties
        (
            IOobject
            (
                "thermophysicalProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );
        thermoType = word(properties.subDict("thermoType").lookup("type"));
    }
