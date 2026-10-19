# 02. Isothermal Ammonia Decomposition

This tutorial introduces heterogeneous ammonia (NH3) decomposition over a nickel catalyst in a porous packed bed under isothermal operation.

The case is simulated using `porousPBRFoam`. The primary objective is to examine the coupling between porous-media flow, species transport, and detailed surface chemistry for NH3 decomposition into nitrogen and hydrogen while the energy equation is disabled.

## Running the case

Load the intended OpenFOAM 8 environment and compile this package as described in the repository installation instructions so that `porousPBRFoam` and the required surface-chemistry libraries are available.

From the repository root:

```bash
cd tutorials/02_Isothermal_NH3Decomposition
./Allrun
```

The `Allrun` script executes the following commands sequentially:

```text
blockMesh
    ↓
porousPBRFoam
```

The script generates the computational mesh and runs the solver selected by `application` in `system/controlDict`. Output is recorded in `log.blockMesh` and `log.porousPBRFoam`. The supplied run controls specify an end time of `5 s` and a fixed time step of `5e-5 s`. The case uses `startFrom latestTime`, so existing time directories can affect the restart point.

To remove the generated mesh and simulation results:

```bash
./Allclean
```

## Files to examine

The following files are particularly useful for understanding the case setup:

- `0/NI(s)` and `0/Ysdefault` — initial surface-species fields
- `constant/porousMediaProperties` — solver mode, isothermal switch, bed properties, and pressure-drop and diffusivity closures
- `constant/surfaceChemistryProperties` — surface-to-volume ratio, surface-chemistry method, and ODE settings
- `constant/reactions.surf` — detailed heterogeneous reaction mechanism
- `constant/catalystThermophysicalProperties` and `constant/thermo.catalyst.solid` — catalyst thermophysical properties
