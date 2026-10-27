const volScalarField& Intr_psiF = thermo.psi();
//const volScalarField& T = thermo.T();
volScalarField& Intr_TF = thermo.T();  //for catalyst
const label inertIndex(composition.species()[inertSpecie]);
