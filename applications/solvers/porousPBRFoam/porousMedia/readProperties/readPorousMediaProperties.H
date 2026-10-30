// Common porous-medium geometry and solid-phase properties.
validatePorousCommonInputs(porousMediaDict);

scalar por(readScalar(porousMediaDict.lookup("Porosity")));
scalar dp(readScalar(porousMediaDict.lookup("ParticleDiameter")));
scalar Cps(readScalar(porousMediaDict.lookup("SolidHeatCapacity")));
scalar rhos(readScalar(porousMediaDict.lookup("SolidDensity")));
// Loaded by the surface-chemistry setup in PMSCHEM mode.
scalar fs2v = 0.0;
