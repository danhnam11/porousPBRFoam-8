# porousPBRFoam_integrate

Rebuilt on cpl2-01 from the server's original porousPBRFoam and porousRhoPBRFoam.
Standard OpenFOAM layout: one main .C file, equation/control .H files, and
Make/files plus Make/options. The solver does not call external solver binaries.

## Build

Load the SAME OpenFOAM environment used for the original solvers (of8x2).
Ensure LIB_PBR8_SRC points to this repository's src directory, then run wmake.
No Windows/WSL object files, generated dependency files, or lnInclude links are supplied.
Do not change compiler, MPI, or label size relative to the existing libraries.
The repository's root Allwmake is unchanged.

## Run

Set application porousPBRFoam_integrate in system/controlDict.
The type in constant/thermophysicalProperties selects the numerical path:

* hePsiThermo: the original porousPBRFoam path, PM or PMSCHEM.
* heRhoThermo: the original porousRhoPBRFoam path, PM only.

rhoConst is an equation-of-state choice, not a property of every rho backend.
The rho path retains the original dynamic mesh and tensor-conductivity behavior.
PMSCHEM with heRhoThermo is rejected explicitly rather than silently ignored.
Use `mpirun -np 4 porousPBRFoam_integrate -parallel` after decomposition.
The standalone OpenFOAM postProcess utility should be used for postprocessing.

## Averaging names

Extr_UF: extrinsic fluid velocity (superficial velocity).
Intr_rhoF, Intr_pF, Intr_TF, Intr_YF, Intr_YiF, Intr_heF: intrinsic fluid fields.
Intr_muF, Intr_CpF, Intr_kappaF: fluid properties at the macroscopic state.
Intr_rhoS, Intr_CpS, Intr_kappaS: solid-phase properties.
Intr_psiF and Intr_dpdtF: compressibility and pressure time derivative.
Intr_heGasF: saved fluid energy field.
F = fluid; S = solid. Under LTE, Intr_TF is the common temperature.

Only C++ identifiers are changed. Existing registry names, disk field names,
species names, member APIs, and dictionary keys remain unchanged. U and rho are
references used by standard OpenFOAM headers, not additional field copies.
Effective coefficients, face mass flow rate phi, source terms and geometry
retain their names. Configuration scalars rhos/Cps/kappas retain their names.

Original-source-hashes.json records the original source/build files and
allows confirmation that the two original solvers were left unchanged.
Numerical validation is left to the user, as requested.

## Server build result

Built successfully on cpl2-01 with OpenFOAM-8.x.2, GCC 4.8.5,
double precision and 64-bit labels (linux64Gcc48DPInt64Opt).
The original 34 source/build files were checked by SHA-256 and are unchanged.
The dependency list was regenerated with the explicit header paths and wmake
completed without missing-header diagnostics or restart loops.
No numerical validation cases were executed.

psiTransportCompatibility.H is a local compatibility wrapper for the three
Rho-named helper aliases in the server's psi transport header. It does not
change the actual transport classes or any original library source.
Continuity-error initialization is shared before selecting the backend,
because the standard initialization header has an include guard.

All source/build-input modification dates match the already-working source
date baseline. This avoids the makefile restart loop caused by future-dated
transferred files on the cluster's shared filesystem.

## Common-code refactor

The main file contains common startup, readThermoType.H, and the selection of
psi/solvePsi.H or rho/solveRho.H. Rho-only configuration validation is in
rho/checkSolverMode.H.

common/ holds shared porous dictionary/property input, thermal-closure inputs
and coefficients, momentum resistance creation, porous/transport/species/energy
field creation, momentum-resistance updates, thermal-state inputs, and effective
species-diffusivity updates. The two createFieldsPorousMedia.H and
updateForPorousMedia.H files now include these shared fragments.

Backend-specific thermal selection, the Standard thermal-transport expression,
Wakao model-name/error behavior, the rho conductivity tensor and Pe field,
surface chemistry, pressure correction and time loops are preserved. No model
class hierarchy or runtime-selection-table conversion is introduced.

The expanded C++ token streams of both property-initialization and closure-update
headers were checked against their pre-refactor versions. Only comments,
formatting and redundant parentheses in sqrt(POR+SMALL) are normalized.
Numerical case validation remains with the user.
