# porousPBRFoam

OpenFOAM 8 solver for volume-averaged porous-media flow, species transport,
heat transfer, and detailed heterogeneous surface chemistry.

## Build and run

Load the intended OpenFOAM 8 environment and configure the package environment,
including `LIB_PBR8_SRC`. From the repository root, build the libraries and
applications using `./Allwmake`. If the required libraries are already built,
rebuild this solver using:

```sh
wmake applications/solvers/porousPBRFoam
```

Set `application porousPBRFoam;` in the case's `system/controlDict`.
Run `porousPBRFoam` from the case directory. For a decomposed case, use
`mpirun -np 4 porousPBRFoam -parallel` with the desired process count.
Use compatible compiler, MPI, precision, and label-size settings for the solver
and its libraries. Run numerical validation separately after rebuilding.

## Thermodynamics and operating modes

The `type` entry in `constant/thermophysicalProperties/thermoType` selects
`psi/solvePsi.H` (`hePsiThermo`) or `rho/solveRho.H` (`heRhoThermo`).
The equation of state determines the density model; `heRhoThermo` alone does
not imply constant density. The constant-density tutorial uses `rhoConst`.

In `constant/porousMediaProperties`:

- `solverMode PM` disables heterogeneous surface-chemistry coupling.
- `solverMode PMSCHEM` enables heterogeneous surface-chemistry coupling and is
  supported by the `hePsiThermo` path. The `heRhoThermo` path accepts `PM` only.
- `isothermal off` (default) solves the macroscopic energy equation.
  `isothermal on` retains the prescribed temperature and skips the thermal
  conductivity closure input and update.

Configure gas-phase chemistry separately through its own dictionaries.

## Source layout

- `porousPBRFoam.C`: startup and thermodynamic-path selection.
- `readThermoType.H`: thermodynamic-type input.
- `psi/`, `rho/`: thermodynamic-path-specific setup, pressure correction,
  chemistry handling, and solution loops.
- `porousMedia/readProperties/`: common properties and one reader per closure
  quantity: `readMomentumResistance.H`, `readSpeciesDiffusivity.H`, and
  `readEffectiveConductivity.H`. The corresponding calculation files are in
  `updateProperties/`. Thermal model selection and all thermal coefficients are
  read together in `readEffectiveConductivity.H`.
- `porousMedia/porousClosureInputValidation.H`: solver-side checks for common
  properties and momentum, diffusivity, and thermal closure input. Wall-condition
  input checks are local to its boundary-condition source file.
- `porousMedia/createFields/`: geometry, resistance, and transport fields.
- `porousMedia/updateProperties/`: resistance and effective-transport updates.
- `porousMedia/equations/`: shared macroscopic conservation equations.
- `Make/files`, `Make/options`: solver build inputs.

## Averaging notation and effective coefficients

`Extr_UF` denotes the extrinsic fluid-phase velocity (superficial velocity).
The `Intr_` prefix denotes intrinsic phase averages; `F` and `S` identify the
fluid and solid phases. Examples are `Intr_rhoF`, `Intr_TF`, `Intr_YiF`, and
`Intr_kappaS`. These C++ identifiers do not rename the existing disk fields
such as `U`, `p`, `T`, or the species fields.

For the Wakao thermal-conductivity closure:

| Manuscript quantity | Code field | Meaning |
| --- | --- | --- |
| k_eff^0 | `kappaEffStag` | Stagnant effective thermal conductivity |
| k_mix | `kappaEffMix` | Flow-induced thermal-dispersion contribution |
| k_eff | `kappaEff` | Sum of the stagnant and mixing contributions |

`StagnantThermalConductivityClosureModel` selects `Kunii`, `Krupiczka`, or
`ZS`. The code evaluates `kappaEffMix = Cmix*Intr_kappaF*Re*Pr` and
`kappaEff = kappaEffStag + kappaEffMix`.

Thermal correlations use the separately prescribed scalar
`SuperficialVelocity` to calculate the particle Reynolds number with the
current fluid properties. This scalar is not automatically read from an inlet
boundary condition. The BruggemanWilhelm diffusivity closure instead uses the
magnitude of the calculated `Extr_UF` field.

## Closure-model inputs

All keywords below belong in `constant/porousMediaProperties` and are
case-sensitive. `Porosity`, `ParticleDiameter`, `SolidDensity`, and
`SolidHeatCapacity` are read for every model selection. The following table
lists additional inputs.

| Quantity | Model | Additional required keywords |
| --- | --- | --- |
| Momentum-exchange source | `None` | None |
| Momentum-exchange source | `Darcy` | `DarcyCoeff` |
| Momentum-exchange source | `Forchheimer` | `DarcyCoeff`, `ForchCoeff` |
| Momentum-exchange source | `Ergun` | None; resistance is calculated from porosity and particle diameter |
| Effective species diffusivity | `BruggemanWilhelm` | None |
| Effective thermal conductivity | `Const` | `EffectiveThermalConductivity` |
| Effective thermal conductivity | `ConstAnisotropic` | `EffectiveRadialThermalConductivity`, `EffectiveAxialThermalConductivity`, `AxialDirection` |
| Effective thermal conductivity | `Standard` | `SolidThermalConductivity` |
| Effective thermal conductivity | `Bunnell` | `SuperficialVelocity` |
| Effective thermal conductivity | `Yagi` | `SolidThermalConductivity`, `SuperficialVelocity`, `YagiAlphaBeta`, `YagiPhi` |
| Effective thermal conductivity | `Wakao` | `SolidThermalConductivity`, `SuperficialVelocity`, `StagnantThermalConductivityClosureModel` |
| Stagnant conductivity within `Wakao` | `Kunii` | `KuniiPhi`, `KuniiGamma` |
| Stagnant conductivity within `Wakao` | `Krupiczka`, `ZS` | None beyond the Wakao inputs |
| Wall-to-bed heat-transfer coefficient | `Const` | `SolidThermalConductivity`, `heatTransferCoefficient` |
| Wall-to-bed heat-transfer coefficient | `Yagi`, `Li` | `SolidThermalConductivity`, `SuperficialVelocity` |
| Wall-to-bed heat-transfer coefficient | `Bey` | `SolidThermalConductivity`, `SuperficialVelocity`, `ReactorDiameter` |

Select the closures using `pressureDropClosureModel`,
`DiffusivityClosureModel`, `ThermalConductivityClosureModel`, and
`WallHeatTransferClosureModel`, respectively. Thermal-conductivity inputs are
read when the energy equation is enabled. Wall-closure inputs apply to the
`porousMediaWallHeatTransfer` boundary condition in `coefficient` mode; its
boundary-field settings are also required. Particle diameter is already a
common input.

`WakaoMixCoeff` is optional and defaults to 0.1. `AxialDirection` must be a
nonzero vector and is normalized internally. The Bey wall correlation uses
the stagnant conductivity supplied by the Wakao closure; select Wakao when
using Bey. Incompatible selections terminate with a model-pairing error.

## Input validation and fatal-error diagnostics

The solver validates common inputs and selected momentum, species-diffusivity,
and thermal-conductivity models before reading their individual coefficients.
The `porousMediaWallHeatTransfer` boundary condition validates its selected
wall closure before reading wall-model inputs in `coefficient` mode.

For missing model inputs, the fatal message prints:

- The selected model or common-property input group.
- The complete required keyword list for that input group.
- Corresponding parameter symbols and units in plain text, where applicable.
- All missing keywords in that group and the dictionary to edit.

Unknown model selections print the available choices. Required scalar values
are checked for finiteness and their permitted sign. Porosity must satisfy
`0 < Porosity < 1`; particle diameter, solid density, solid heat capacity,
and conductivities must be positive. Resistance coefficients, reference
velocity, mixing coefficients, and a prescribed wall coefficient may be zero
but must not be negative. The anisotropic axial direction must be finite and
nonzero and is normalized by the reader.

For Wakao, the stagnant-conductivity model is mandatory. Kunii additionally
requires `KuniiPhi` and `KuniiGamma`. The optional `WakaoMixCoeff` retains its
default of 0.1. Bey requires `ThermalConductivityClosureModel Wakao` and a
reactor diameter greater than the particle diameter; incompatible selections
terminate instead of using a zero stagnant-conductivity contribution.

Validation stops at the first failing group (for example, common inputs before
model-specific inputs). Incorrect value types or malformed OpenFOAM syntax
remain subject to the standard dictionary parser's error messages. These
checks do not establish whether a correlation is applicable to a particular
experiment or operating regime.

Solver-side checks are implemented in
`applications/solvers/porousPBRFoam/porousMedia/porousClosureInputValidation.H`.
Wall-model checks are local to
`src/catalyticSurfaceBoundaries/porousMediaWallHeatTransfer/porousMediaWallHeatTransferFvPatchScalarField.C`.
The boundary-condition library does not depend on solver headers. Rebuild the
solver after changing its checks, and the boundary-condition library after
changing wall checks. The repository-root `./Allwmake` builds both in order.

The validation helper was compiled against OpenFOAM 8 and exercised with 41
valid, missing-input, invalid-value, and incompatible-model cases. This check
does not replace a full solver build or numerical case validation. After the
reader reorganization, 36 combinations (including isothermal operation) produced
the same initialized values as the original reader headers; zero axial-vector
rejection and the 41 validation cases were also checked.

### Integrated build and startup checks (2026-10-30)

The repository-root `./Allwmake -j4` completed successfully with OpenFOAM 8
(`linux64GccDPInt32Opt`). Temporary copies of the wall-heated non-reacting
(`heRhoThermo`) and isothermal NH3 decomposition (`hePsiThermo`, `PMSCHEM`)
tutorials each completed three time steps using the rebuilt solver. The
reacting case also completed surface-mechanism conversion.

The rebuilt solver was exercised with six deliberately invalid inputs:
missing Forchheimer coefficient, missing Yagi coefficients, missing Wakao
stagnant-model selection, missing Kunii coefficients, missing constant wall
coefficient, and an incompatible Bey/non-Wakao selection. Each terminated
with the expected model-specific fatal diagnostic.

These are build and short execution checks. They do not establish steady-state
convergence or reproduce the manuscript validation results. Original tutorial
results were not overwritten.

### Fatal-message formatting

Missing-input diagnostics use a single table with `Keyword`, `Symbol`, `Unit`,
and `Status` columns. Each column is left-aligned with a fixed width calculated
from its longest entry, and status is shown as `provided` or `MISSING`. Model,
selection keyword, dictionary path, and the action to take surround the table.
Unsupported models and invalid scalar values use the same labelled layout.
The standard OpenFOAM fatal banner and source-location report are retained;
no ANSI colour sequences are added. The input-validation conditions are
unchanged. After this formatting-only revision, the boundary-condition library
and solver were rebuilt successfully; no cases were executed for this revision.
