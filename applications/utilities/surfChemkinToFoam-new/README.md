# surfChemkinToFoam-new

Independent OpenFOAM Foundation 8 surface-mechanism converter. The original
`surfChemkinToFoam` directory, executable and tutorial scripts are not replaced.
The executable and production sources use C++11 and require no Python runtime,
generated lexer, or external chemistry package. Python is used only by the
optional regression checks.

## Build

Source your OpenFOAM 8 environment and set `LIB_PBR8_SRC` to this repository's
`src` directory, as for the existing utility. Then run `wmake` in this directory.
The resulting executable is `$FOAM_USER_APPBIN/surfChemkinToFoam-new`.
This directory is intentionally not added to the repository's existing Allwmake.

## Existing invocation

The original five positional arguments are unchanged:

```sh
surfChemkinToFoam-new \
    surf.inp surfTherm.dat surfaceThermoDict reactions.surf thermo.solid
```

Original constant-property defaults are supported unchanged:

```text
surfaceThermoDefaults
{
    rho 8900;
    mu 0;
    Pr 1;
    Cp 440;
    Cv 440;
    Hf 0;
}
```

`rho`, `mu` and `Pr` remain required. Cp, Cv and Hf default to zero as in the
original converter. These values are user inputs; zero is not inferred physical
thermochemistry. In particular, ensure formation enthalpies are appropriate
when interpreting catalytic heat release.

## Additional features

### Validate without writing

```sh
surfChemkinToFoam-new -checkOnly \
    surf.inp surfTherm.dat surfaceThermoDict reactions.surf thermo.solid
```

Checks include SITE definitions, positive site density, numeric syntax,
NASA7 coefficient-record structure, missing required species thermo,
mapping collisions, site conservation, modifiers, property values, and
input/output path collisions. Errors include the input line where available.
All data are validated before output files are replaced. Output files are staged
in their destination directories; previous outputs are restored if installation
of either file fails during normal error handling.

This is not a kinetic validation or an element/mass-balance validator. The five
original inputs do not include a separate gas-species registry. Multiple SITE
phases and non-unit site occupancy are rejected rather than flattened.

### Input units

The REACTIONS header accepts CAL/MOLE, KCAL/MOLE, JOULES/MOLE, KJOULES/MOLE,
KELVINS, or EVOLTS, and MOLES or MOLECULES. D/d exponent notation is accepted.
Concentration prefactors retain CHEMKIN centimetre-based length units and are
converted to the existing OpenFOAM kmol/m3 and kmol/m2 convention. SDEN is read
in mol/cm2 and converted to kmol/m2. STICK prefactors are not concentration-scaled.
Coverage epsilon uses the selected activation-energy unit.

If the energy unit is absent, the old KJOULES/MOLE assumption is retained with
a warning; this differs from generic CHEMKIN defaults, so specify units explicitly.
MWOFF is accepted. Other unimplemented header options are rejected.

### Bare-site names

```sh
surfChemkinToFoam-new -siteSpecies 'Ni(s)' \
    surf.inp surfTherm.dat surfaceThermoDict reactions.surf thermo.solid
```

This maps `*` to `Ni(s)` and names such as `H(*)` to `H(s)` consistently in SITE,
reaction equations, COV entries, and thermo species headers. It does not modify
the input files. SITE species must appear on lines after the SITE/SDEN header.

Alternatively, add an explicit map to surfaceThermoDict:

```text
speciesMap
{
    "*"     "Ni(s)";
    "H(*)"  "H(s)";
    "NI(s)" "Ni(s)";
}
```

Mappings are exact, one-pass, and case-sensitive. Explicit mappings take
precedence over -siteSpecies. Targets must be unique for distinct SITE species.
The matching thermo names must be mapped to the same targets.

### Per-species constant properties

Add optional overrides to surfaceThermoDict. Unspecified values inherit defaults.
Use the final mapped species names as keys:

```text
surfaceThermoOverrides
{
    "Ni(s)"
    {
        rho 8900;
        Cp 440;
    }
    "H(s)"
    {
        Cp 123;
        Hf -456;
    }
}
```

Allowed properties: rho, mu, Pr, Cp, Cv, Hf, molWeight. rho, Pr and molWeight
must be positive; mu, Cp and Cv must be nonnegative. Hf may have either sign.
An explicit molWeight is required for labels whose chemical formula cannot be
inferred. As in the original converter, adsorbate molecular weight is inferred
from the species label without adding a substrate atom to every adsorbate.
Bare uppercase metal labels NI(s), PT(s), RH(s), RU(s) are handled as Ni, Pt,
Rh, Ru; the old NI(s) -> N + I molecular-weight error is corrected.

### Summary and output

The converter reports species and reaction counts, STICK/COV/PKC counts, energy
and amount units, site density, mappings, and warnings. Both output object names
use the requested filename. Only SITE species receive thermo property blocks;
unrelated thermo records are not written. Repeated unused records are tolerated
with a warning. A complete thermo record at EOF without END is accepted with a
legacy-compatibility warning; truncated coefficient records are rejected.

## Retained scope and explicit limitations

- Ordinary Arrhenius, STICK, COV and PKC retain the existing runtime reaction types.
- Only explicit irreversible `=>` reactions are supported. Reversible `=`/`<=>`,
  REV, FORD, DUP, BULK and unknown modifiers are not silently approximated.
- Combined STICK+COV or PKC+COV requires an appropriate combined library rate
  model and is rejected here.
- Constant-property thermo output is retained. NASA7 coefficient records are
  read and checked but do not automatically replace Cp or Hf. NASA9 is unsupported.
- `-newFormat` is retained as a deprecated command-line compatibility option;
  NASA7 record reading is automatic, not controlled by this flag.

## Verification

With the OpenFOAM environment loaded:

```sh
python3 testTutorials.py
python3 testFeatures.py
```

The repository root is inferred from this installed source directory. For a
different layout, set POROUS_PBR_REPO. Checks use temporary directories and do
not overwrite tutorial inputs or outputs.

Regression comparison covers all nine current conversion tutorials, 196
reactions: reaction types, equations, Arrhenius and coverage coefficients match
the old utility. Required SITE species and constant properties also match,
except the intentional correction of uppercase NI(s) molecular weight.
The new thermo output omits unused extra property blocks from the input database.

Feature checks cover unit equivalence, molecule units, species mapping,
per-species overrides, validation-only operation, malformed/unsupported inputs,
preservation of existing outputs on validation errors, and reading the generated
dictionaries with foamDictionary. This does not substitute for reactor-simulation
validation or a separate build on the older cluster toolchain.
