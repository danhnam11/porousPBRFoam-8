const volScalarField& psi = thermo.psi();
//const volScalarField& T = thermo.T();
volScalarField& T = thermo.T();  //for catalyst
const label inertIndex(composition.species()[inertSpecie]);
