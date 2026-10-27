# 00. Surface CHEMKIN Conversion

This tutorial introduces the conversion of CHEMKIN-format surface mechanisms into the OpenFOAM dictionaries used by the surface-chemistry library.

The conversion is performed using `surfChemkinToFoam`. The primary objective is to become familiar with the mechanism inputs, surface thermophysical settings, and generated reaction and species dictionaries before setting up a catalytic reactor simulation.

## Running the case

Load the intended OpenFOAM 8 environment and compile this package as described in the repository installation instructions so that `surfChemkinToFoam` is available.

From the repository root, select an example subdirectory:

```bash
cd tutorials/00_surfChemkinToFoam/example
./Allrun
```

The `Allrun` script executes the following commands sequentially:

```text
Allclean
   ↓
surfChemkinToFoam
```

The conversion command is:

```bash
surfChemkinToFoam \
    surf.inp \
    surfTherm.dat \
    surfaceThermoDict \
    reactions.surf \
    thermo.solid
```

Run the script from the selected subdirectory. Its initial cleanup removes the previous `reactions.surf` and `thermo.solid` before regenerating them. To use another mechanism, enter its subdirectory and run its `Allrun` script instead.

To remove the generated dictionaries:

```bash
./Allclean
```

## Files to examine

The following files in each example subdirectory are particularly useful for understanding the conversion:

- `surf.inp` — CHEMKIN-format surface species, site information, and reaction mechanism
- `surfTherm.dat` — thermodynamic data supplied to the converter
- `surfaceThermoDict` — default surface density, transport, and constant heat-capacity settings
- `reactions.surf` — generated OpenFOAM surface reaction dictionary
- `thermo.solid` — generated constant-property surface species dictionary
- `Allrun` — cleanup and conversion commands
