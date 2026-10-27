# 00. Surface CHEMKIN Conversion

This tutorial introduces the conversion of CHEMKIN-format surface mechanisms into the OpenFOAM dictionaries used by the surface-chemistry library.

The conversion is performed using `surfChemkinToFoam`. The primary objective is to become familiar with the mechanism inputs, surface thermophysical settings, and generated reaction and species dictionaries before setting up a catalytic reactor simulation.

## Running the case

Load the intended OpenFOAM 8 environment and compile this package as described in the repository installation instructions so that `surfChemkinToFoam` is available.

From the repository root, select an example subdirectory:

```bash
cd tutorials/00_surfChemkinToFoam/2013_Karakaya_Rh_H2Oxidation_example
./Allchemkin
```

The `Allchemkin` script executes the following commands sequentially:

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

Run the script from the selected subdirectory. Its initial cleanup removes the previous `reactions.surf` and `thermo.solid` before regenerating them. To use another mechanism, enter its subdirectory and run its `Allchemkin` script instead.

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
- `Allchemkin` — cleanup and conversion commands

## Mechanism directories

Directory names use publication year, first author, catalyst and reaction. AKC and PKC retain the original formulation labels. The example directory is an unchanged duplicate of the Karakaya H2-oxidation input set.

- `2011_Appari_Ni_NH3Decomposition`
- `1996_Deutschmann_Pt_CH4Oxidation`
- `2013_Karakaya_Rh_H2Oxidation`
- `2014_Karakaya_Rh_WGS`
- `2011_Maier_Ni_SMR`
- `2021_Schmider_Ni_CO_CO2Methanation`
- `2016_Takahashi_Ru_NH3Decomposition_AKC`
- `2016_Takahashi_Ru_NH3Decomposition_PKC`
- `2013_Karakaya_Rh_H2Oxidation_example`

Karakaya H2 oxidation uses the 2013 journal publication year (DOI: 10.1016/j.ces.2012.11.004), rather than the old 2012 label. The 2021 methanation mechanism is attributed to first author Schmider (DOI: 10.1021/acs.iecr.1c00389). Names identify mechanism provenance; they do not validate modified input files.
