#include "porousMedia/updateProperties/updateMomentumResistance.H"

#include "porousMedia/updateProperties/updateEffectiveConductivity.H"

// Update the effective species diffusivities using the selected closure.
if (DiEffClosureModel_ == "BruggemanWilhelm")
{
#include "porousMedia/updateProperties/updateSpeciesDiffusivity.H"


}
