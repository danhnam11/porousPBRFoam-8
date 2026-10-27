// Reset catalytic heat release for PM and refresh it for PMSCHEM.
QdotSurface = dimensionedScalar("zero", QdotSurface.dimensions(), 0.0);

if (surfaceReaction.valid())
{
    surfaceReaction->solveSurfaceReactions();
    QdotSurface = FS2V*surfaceReaction->QdotFromSCHEM();
}
