# 02. Isothermal Ammonia Decomposition

This tutorial introduces heterogeneous ammonia (NH3) decomposition over a nickel catalyst in a porous packed bed under isothermal operation.

The case uses `porousPBRFoam` to demonstrate the coupling between porous-media flow, species transport, and detailed surface chemistry for NH3 decomposition into nitrogen and hydrogen. The `solverMode` is set to `PMSCHEM` and `isothermal` is set to `on` in `constant/porousMediaProperties`, so the energy equation is disabled and the prescribed temperature field is retained.

## Running the case

Load the intended OpenFOAM 8 environment and compile this package as described in the repository installation instructions so that `porousPBRFoam`, `surfChemkinToFoam`, and the required surface-chemistry libraries are available.

From the repository root:

```bash
cd tutorials/02_Isothermal_NH3Decomposition
./Allrun
```

The `Allrun` script performs the following steps in sequence:

1. Runs `chemkin/Allchemkin` to convert the CHEMKIN surface mechanism and thermodynamic data using `surfChemkinToFoam`, generating `constant/reactions.surf` and `constant/thermo.solid`.
2. Runs `blockMesh` to generate the computational mesh.
3. Runs the solver selected by `application` in `system/controlDict`, which is `porousPBRFoam` in this tutorial.

Mesh-generation and solver output are recorded in `log.blockMesh` and `log.porousPBRFoam`, respectively. Mechanism-conversion output is printed to the terminal.

## Cleaning the case

To remove the generated mesh, simulation results, and converted surface-chemistry files:

```bash
./Allclean
```

The converted files are regenerated the next time `./Allrun` is executed.

## Files to examine

- `0/Ni(s)` and `0/Ysdefault`: initial surface-species fields.
- `constant/porousMediaProperties`: solver mode, isothermal switch, bed properties, and momentum-resistance and diffusivity closure models.
- `constant/surfaceChemistryProperties`: surface-to-volume ratio, surface-chemistry method, and ODE integration settings.
- `chemkin/Allchemkin`: surface-mechanism conversion script.
- `chemkin/surf.inp`, `chemkin/surfTherm.dat`, and `chemkin/surfaceThermoDict`: inputs to `surfChemkinToFoam`.
- `constant/reactions.surf`: generated heterogeneous reaction mechanism.
- `constant/catalystThermophysicalProperties` and `constant/thermo.solid`: surface-species thermophysical configuration and generated property data.
