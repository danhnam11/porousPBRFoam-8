# porousPBRFoam-8

## General Information
`porousPBRFoam-8` is an OpenFOAM 8 solver for volume-averaged porous-media simulations of catalytic packed-bed reactors. The solver couples macroscopic fluid-flow, species-transport, and energy equations with detailed gas-phase and heterogeneous surface chemistry.

The solver provides runtime-selectable closure models for:

- Pressure drop through packed beds
- Effective thermal conductivity
- Effective species diffusivity
- Heat transfer between the reactor wall and the porous bed

Heterogeneous reaction rates evaluated per unit catalytic surface area are converted into volumetric species and heat source terms using the specified catalytic surface-to-volume ratio.

The surface-chemistry capability is based on the [`surfaceChemistryFoam-8`](https://github.com/danhnam11/surfaceChemistryFoam-8) library.

## Installation
- The complete installation of the OpenFOAM 8.0 framework in a Linux operating system is required before installing this package, as it is designed for the Linux-based OpenFOAM 8.0 version. 
- Prepare a directory on your system, for example, _yourDirectory_:

		mkdir ~/OpenFOAM/yourDirectory/
		cd ~/OpenFOAM/yourDirectory/	
- Download source files using git: 

		git clone https://github.com/danhnam11/porousPBRFoam-8.git

- Specify the path of the _src_ directory of this package to an environment variable named _LIB_PBR8_SRC_. Suppose the _porousPBRFoam-8_ have been downloaded into _yourDirectory_. Then the following commands should be executed to specify the path of the _src_:

		echo "export LIB_PBR8_SRC=~/OpenFOAM/yourDirectory/porousPBRFoam-8/src/" >> ~/.bashrc
		source ~/.bashrc

- To compile the necessary libraries and solver, go to _porousPBRFoam-8_ directory and run the _Allwmake_ script (it may take hours to be finished):

		cd ~/OpenFOAM/yourDirectory/porousPBRFoam-8/
		./Allwmake

- After successful compilation, the following libraries are saved at _$FOAM_USER_LIBBIN_ :

		libcatSpecie.so
		libthermophysicalProperties.so
		libfluidThermophysicalModels.so
		libcatFluidThermophysicalModels.so		
		libreactionThermophysicalModels.so
		libcatReactionThermophysicalModels.so		
		libchemistryModel.so		
		libsurfaceChemistryModel.so		
		libsolidThermo.so
		libSLGThermo.so
		libfluidThermoMomentumTransportModels.so
		libthermophysicalTransportModels.so
		libpsiReactionThermophysicalTransportModels.so
		librhoReactionThermophysicalTransportModels.so
		libradiationModels.so
		libcombustionModels.so
		libfvOptions.so
		libspecieTransfer.so		
		libcatalyticSurfaceFvPatchFields.so
		libsurfaceFilmModels.so
		libsurfaceFilmDerivedFvPatchFields.so
		libthermalBaffleModels.so
		libfieldFunctionObjects.so
		libforces.so
		liblagrangianFunctionObjects.so
		libsolverFunctionObjects.so
		libutilityFunctionObjects.so
		
- and the following executable programs are saved at _$FOAM_USER_APPBIN_ :

		catalystFoam
		porousPBRFoam
		surfChemkinToFoam
		
		DTLreactingFoam
		DTMchemkinToFoam
		FTMchemkinToFoam
		
		reactingFoam
		chemkinToFoam		


- These newly compiled libraries, solvers, and utilities are now ready for use.
- It is important to note that if a different solver (i.e., program) relies on any of the aforementioned compiled libraries, its corresponding _options_ file, located in the _Make_ directory, must be updated accordingly. The solver should then be recompiled to prevent potential conflicts, such as segmentation faults. For reference, the _Make_ directory of the _reactingFoam_ solver included in this package provides a convenient example. Although this version of _reactingFoam_ is nearly identical to the original, its _options_ file has been modified to link against compiled libraries in this package. As a result, it can be used seamlessly in place of the original version without issue.

- To remove all compiled libraries and solvers, go to _porousPBRFoam-8_ directory and run the _Allwclean_ script:

		cd ~/OpenFOAM/yourDirectory/porousPBRFoam-8/
		./Allwclean

## Using this package 
Upon completing the compilation process, the _porousPBRFoam_ solver can be utilized by simply typing its name in the terminal, _porousPBRFoam_. All important instructions for using surface chemistry library in a new developed OF-based solver and case setting are provided in _documentations_ directory.

## Tutorials
A test cases in catalytic processes is available in the _tutorials_ directory.

	cd ~/OpenFOAM/yourDirectory/porousPBRFoam-8/tutorials/

## Authors 
This package was developed at the Clean Combustion & Energy Research Lab., Dept. of Mech. Engineering, Ulsan National Institute of Science and Technology (UNIST), Korea (Prof. C.S. Yoo: https://csyoo.unist.ac.kr/). If you publish results obtained by using this package, please cite our paper as follows:
- J. H. Lee, D. N. Nguyen, H. W. Seo, H. J. Ahn, C. S. Yoo, Volume-averaged porous-media modeling of catalytic packed-bed reactors with detailed surface chemistry: Implementation and validation of porousPBRFoam, Chemical Engineering Journal (2026)(submitted).

Contact:
- danhnam11@gmail.com or jhlee28@unist.ac.kr

## Reference
