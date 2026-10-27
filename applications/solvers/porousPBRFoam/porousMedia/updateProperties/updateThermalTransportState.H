Re = (DP*USUP)*Intr_rhoF/Intr_muF;
Pr = Intr_CpF*Intr_muF/Intr_kappaF;

// Initialize the effective thermal-transport fields using fluid properties.
kappaEff = Intr_kappaF;
pAlphaEff = thermophysicalTransport->alphaEff();

// Update the effective thermal conductivity using the selected closure.
