# 01. Wall-Heated Packed-Bed Heat Transfer

This tutorial introduces non-reacting heat transfer in a wall-heated porous packed bed.

The case is simulated using `porousRhoPBRFoam`, a solver designed for plug-flow modeling under the assumption of constant fluid density, $\rho_f = \mathrm{const}$. The primary objective is to examine the porous-media thermal-conductivity and wall heat-transfer settings without the additional complexity of chemical reactions.

## Running the case

Load the intended OpenFOAM 8 environment and compile this package as described in the repository installation instructions so that `porousRhoPBRFoam` and the required libraries are available.

From the repository root:

```bash
cd tutorials/01_WallHeated_HeatTransfer
./Allrun
```

The `Allrun` script executes the following commands sequentially:

```text
blockMesh
    ↓
porousRhoPBRFoam
```

The script generates the computational mesh and runs the solver selected by `application` in `system/controlDict`. Output is recorded in `log.blockMesh` and `log.porousRhoPBRFoam`. The case uses `startFrom latestTime`, so existing time directories can affect the restart point.

To remove the generated mesh and simulation results:

```bash
./Allclean
```

## Files to examine

The following files are particularly useful for understanding the case setup:

- `constant/porousMediaProperties` — bed geometry, solid properties, and closure-model selections
