//QdotSurface = (POR+SMALL)*FS2V*catReaction->QdotFromSurf();
//Qdot = (POR+SMALL)*reaction->Qdot();

QdotSurface = FS2V*surfaceReaction->QdotFromSCHEM();
Qdot = reaction->Qdot();
