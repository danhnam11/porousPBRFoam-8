// The server's psi transport header uses Rho names for three helper aliases.
// Give those aliases distinct names in this translation unit only. The actual
// transport class and original library files remain unchanged.
#define laminarRhoReactionThermophysicalTransportModel laminarPsiReactionThermophysicalTransportModel
#define RASRhoReactionThermophysicalTransportModel RASPsiReactionThermophysicalTransportModel
#define LESRhoReactionThermophysicalTransportModel LESPsiReactionThermophysicalTransportModel
#include "psiReactionThermophysicalTransportModel.H"
#undef laminarRhoReactionThermophysicalTransportModel
#undef RASRhoReactionThermophysicalTransportModel
#undef LESRhoReactionThermophysicalTransportModel
