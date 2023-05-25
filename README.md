# A curated list of repositories related to fluid dynamics.
## 1. Learning Materials

- [barbagroup/CFDPython](https://github.com/barbagroup/CFDPython): A sequence of Jupyter notebooks featuring the "12 Steps to Navier-Stokes".
- [gpeyre/numerical-tours](https://github.com/gpeyre/numerical-tours): Numerical Tours of Signal Processing and other materials.
- [jfavre/Visualization-training](https://github.com/jfavre/Visualization-training): The material used for the Scientific Visualization course, organized online by the Swiss National Supercomputing Centre (CSCS) on May 17-18, 2021.
- [surajp92/CFD_Julia](https://github.com/surajp92/CFD_Julia): This repository contains fundamental codes related to CFD that can be included in any graduate level CFD coursework.
- [cdegroot/cfdcourse](https://bitbucket.org/cdegroot/cfdcourse): A course on Computational Fluid Dynamics using Jupyter Notebooks and Python.
- [Steve Brunton/Fluid Dynamics](https://www.youtube.com/playlist?list=PLMrJAkhIeNNQWO3ESiccZmPssvUDFHL4M): Prof. Brunton's lecture series on Fluid dynamics.
- [Barry Belmont/NSF Fluid Mechanics Series](https://www.youtube.com/playlist?list=PL0EC6527BE871ABA3): A collection of NFS Fluid Mechanics lecture series from mid 20th century.
- [TME226 Mechanics of fluids](http://www.tfd.chalmers.se/~lada/MoF/): Fluid Mechanics course by Prof. Lars Davidson.
- [CH EN 6355 – Computational Fluid Dynamics](http://www.tonysaad.net/ucfd/): Computational Fluid Dynamics course by Prof. Tony Saad.
- [CSE-lab Teaching](https://www.cse-lab.ethz.ch/teaching/): Course about high performance computing by CSE-lab directed by Prof. Petros Koumoutsakos.
- [Deep Bayes Course](http://deepbayes.ru/2019/#): Summer School on Deep Learning and Bayesian Methods.
- [An Introduction to Statistical Learning](https://www.statlearning.com/online-course): A free online companion course to the Second Edition of An Introduction to Statistical Learning.
- [CS106B](http://web.stanford.edu/class/archive/cs/cs106b/cs106b.1184/lectures.shtml): Programming Abstractions course.
- [Programming in C/C++ Efficiently](https://github.com/ShiqiYu/CPP): Course 'CS205 C/C++ Program Design' in 2021 Fall at Southern University of Science and Technology.

### Books

- Stephen B. Pope. [Turbulent Flows](https://pope.mae.cornell.edu/TurbulentFlows/solutions/index.html)
- C. Hirsch. [Numerical Computation of Internal and External Flows](https://waxworksmath.com/Authors/G_M/Hirsch/hirsch.html)
- Katate Masatsuka. [I do like CFD](http://ossanworld.com/cfdbooks/index.html)
- Doyub Kim. [Fluid Engine Development](https://fluidenginedevelopment.org/)
- Joel H. Ferziger, Milovan Peric, Robert L. Street. [Computational Methods for Fluid Dynamics](http://www.cfd-peric.de/).

### Blogs

- [The Visual Room](http://www.thevisualroom.com/index.html)
- [All About CFD](https://cfdisraelblog.wordpress.com/)
- [A complete path to learn CFD](https://www.jahid-hasan.com/writings/a-complete-learning-path-for-cfd/)
- [CFD WITH A MISSION](https://caefn.com/)

### Researchers

- [Prof. Lars Davidsson](http://www.tfd.chalmers.se/~lada/)
- [Prof. Sandip Mazumder](https://www.youtube.com/@sandipmazumder171)
- [Prof. Steven L. Brunton](https://www.eigensteve.com/people)

## 2. Computational fluid dynamics (CFD)

### Solvers

#### Finite element methods (FEM)

- [JuliaFEM/JuliaFEM.jl](https://github.com/JuliaFEM/JuliaFEM.jl): The JuliaFEM software library is a framework that allows for the distributed processing of large Finite Element Models across clusters of computers using simple programming models. It is designed to scale up from single servers to thousands of machines, each offering local computation and storage.
- [FEniCS/dolfinx](https://github.com/FEniCS/dolfinx): Next generation FEniCS problem solving environment
- [deal.II](https://dealii.org/): An open source finite element library
- [firedrakeproject/firedrake](https://github.com/firedrakeproject/firedrake): Firedrake is an automated system for the portable solution of partial differential equations using the finite element method (FEM)
- [KratosMultiphysics/Kratos](https://github.com/KratosMultiphysics/Kratos):  Kratos Multiphysics (A.K.A Kratos) is a framework for building parallel multi-disciplinary simulation software. Modularity, extensibility and HPC are the main objectives. Kratos has BSD license and is written in C++ with extensive Python interface.
- [feelpp/feelpp](https://github.com/feelpp/feelpp): Finite Element Embedded Language and Library in C++
- [FemusPlatform/femus](https://github.com/FemusPlatform/femus): Multigrid Finite Element Code FEMuS
- [Netgen/NGSolve](https://ngsolve.org/): High performance multiphysics finite element software
- [ngsxfem/ngsxfem](https://github.com/ngsxfem/ngsxfem): Add-On to NGSolve for unfitted finite element discretizations
- [PETSc-FEM](https://cimec.org.ar/foswiki/Main/Cimec/PETScFEM): A general purpose, parallel, multi-physics FEM (Finite Element Method) program for CFD (Computational Fluid Dynamics) applications based on PETSc
- [Firedrake ](https://www.firedrakeproject.org/): Firedrake is an automated system for the solution of partial differential equations using the finite element method
- [FEATFLOW ](http://www.mathematik.tu-dortmund.de/~featflow/en/index.html): High performance finite elements
- [ryujin](https://github.com/conservation-laws/ryujin): High-performance second-order collocation-type finite-element scheme for solving the compressible Navier-Stokes and Euler equations of gas dynamics on unstructured meshes
- [gascoigne](https://kosinus.math.uni-magdeburg.de/gascoigne): a flexible finite element toolkit

#### Finite volume methods (FVM)

- [OpenFOAM/OpenFOAM-dev](https://github.com/OpenFOAM/OpenFOAM-dev): OpenFOAM is a free, open source computational fluid dynamics (CFD) software package released by the OpenFOAM Foundation.
- [su2code/SU2](https://github.com/su2code/SU2): SU2: An Open-Source Suite for Multiphysics Simulation and Design
- [cselab/aphros](https://github.com/cselab/aphros): Finite volume solver for incompressible multiphase flows with surface tension
- [code-saturne/code_saturne](https://github.com/code-saturne/code_saturne):  code_saturne public mirror
- [FiPy](https://www.ctcms.nist.gov/fipy/): A Finite Volume PDE Solver Using Python
- [FLUBIO ](https://gitlab.com/alie89/flubio-code-fvm): An unstructured, parallel, finite-volume based Navier-Stokes and convection-diffusion like equations solver for teaching and research purposes

#### Finite Difference Methods (FDM)

- [pencil-code/pencil-code](https://github.com/pencil-code/pencil-code): A high-order finite-difference code for compressible hydrodynamic flows with magnetic fields and particles
- [roshansamuel/saras](https://github.com/roshansamuel/saras): An MPI parallelized Navier-Stokes equation solver written in C++. It uses the finite-difference method for calculating spatial derivatives and parallelized geometric multi-grid method for solving the pressure Poisson equation
- [OpenSBLI](https://opensbli.github.io/): A framework for the automated derivation of finite difference solvers from high-level problem descriptions
- [xcompact3d/Incompact3d](https://github.com/xcompact3d/Incompact3d): High-order finite-difference flow solvers dedicated to the study of turbulent flows

#### Spectral methods

- [DedalusProject/dedalus](https://github.com/DedalusProject/dedalus):  A flexible framework for solving PDEs with modern spectral methods.
- [FourierFlows/FourierFlows.jl](https://github.com/FourierFlows/FourierFlows.jl): Tools for building fast, hackable, pseudospectral partial differential equation solvers on periodic domains
- [Nek5000/Nek5000](https://github.com/Nek5000/Nek5000): NEK5000 is an spectral element CFD code developed at the Mathematics and Computer Science Division of Argonne National Laboratory.
- [spectralDNS/shenfun](https://github.com/spectralDNS/shenfun): High performance computational platform in Python for the spectral Galerkin method
- [Semtex](https://users.monash.edu.au/~bburn/semtex.html): Semtex is a classical quadrilateral spectral element incompressible direct numerical simulation code that uses the standard nodal GLL basis functions to provide two-dimensional solutions and (optionally) Fourier expansions in a homogeneous direction to provide three-dimensional solutions
- [Nek5000/nekRS](https://github.com/Nek5000/nekRS): Navier Stokes solver based on the spectral element method targeting classical processors and hardware accelerators like GPUs
- [neko](https://github.com/ExtremeFLOW/neko): Neko is a portable framework for high-order spectral element flow simulations
- [channelflow](https://github.com/epfl-ecps/channelflow): Channelflow is a software system for numerical analysis of the incompressible fluid flow in channel geometries, written in C++ and MPI-parallelized

#### Smoothed Particle Hydrodynamics (SPH) Methods

- [Xiangyu-Hu/SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys):  It provides C++ APIs for physical accurate simulation and aims to model coupled industrial dynamic systems including fluid, solid, multi-body dynamics and beyond with SPH (smoothed particle hydrodynamics), a meshless computational method using particle discretization
- [DualSPHysics](https://dual.sphysics.org/): DualSPHysics is based on the Smoothed Particle Hydrodynamics model named SPHysics (www.sphysics.org). The code is developed (GNU Lesser General Public License) to study free-surface flow phenomena where Eulerian methods can be difficult to apply. DualSPHysics is a set of C++ and CUDA codes designed to deal with real-life engineering problems
- [SPlisHSPlasH](https://github.com/InteractiveComputerGraphics/SPlisHSPlasH): SPlisHSPlasH is an open-source library for the physically-based simulation of fluids. The simulation in this library is based on the Smoothed Particle Hydrodynamics (SPH) method which is a popular meshless Lagrangian approach to simulate complex fluid effects
- [pypr/pysph](https://github.com/pypr/pysph): A framework for Smoothed Particle Hydrodynamics in Python
- [SPHysics](https://wiki.manchester.ac.uk/sphysics/index.php/Main_Page): SPH Free-surface Flow Solver

#### Lattice Boltzmann methods (LBM)

- [aromanro/LatticeBoltzmann](https://github.com/aromanro/LatticeBoltzmann): A 2D Lattice Boltzmann program
- [Palabos](https://palabos.unige.ch/): The Palabos library is a framework for general-purpose computational fluid dynamics (CFD), with a kernel based on the lattice Boltzmann (LB) method
- [lbmpy](https://i10git.cs.fau.de/pycodegen/lbmpy/): Run fast fluid simulations based on the lattice Boltzmann method in Python on CPUs and GPUs
- [OpenLB](https://www.openlb.net/): The OpenLB project provides a C++ package for the implementation of lattice Boltzmann methods that is general enough to address a vast range of tansport problems
- [TCLB](https://github.com/CFD-GO/TCLB): Templated MPI+CUDA/CPU Lattice Boltzmann code
- [hemelb](https://github.com/hemelb-codes/hemelb): A high performance parallel lattice-Boltzmann code for large scale fluid flow in complex geometries

#### Immersed Boundary Methods (IBM)
- [FluSI & Wabbit](http://aifit.cfd.tu-berlin.de/wordpress/?page_id=42): Aerodynamics of Insect Flight and Turbulent Flow
- [LESGO](https://lesgo.me.jhu.edu/): LESGO is a parallel pseudo-spectral large-eddy simulation code
- [notus](https://notus-cfd.org/): It is dedicated to the modelisation and simulation of incompressible fluid flows in a massively parallel context. Its numerical framework is the Finite Volume method on Cartesian staggered grids with a methodological focus on interfaces treatment
- [Lethe](https://github.com/lethe-cfd/lethe): Open-source computational fluid dynamics (CFD) software which uses high-order continuous Galerkin formulations to solve the incompressible Navier–Stokes equations
- [OpenIFEM/OpenIFEM](https://github.com/OpenIFEM/OpenIFEM): An implementation of the Immersed Finite Element Method based on deal.II
- [stochasticHydroTools/IBMethod](https://github.com/stochasticHydroTools/IBMethod): Immersed Boundary Method
- [Virtual Flow Simulator](https://zenodo.org/record/4677354): The Virtual Flow Simulator (VFS-Geophysics) code is a curvilinear immersed boundary-based model. It includes DNS, K-omega, and LES turbulence models
- [YunchaoYang/NekIBM](https://github.com/YunchaoYang/NekIBM): A Multiphase flow simulation platform using Direct-forced Immersed Boundary Method based on Spectral element solver Nek5000

#### Discontinuous Galerkin Methods (DG)

- [exadg/exadg](https://github.com/exadg/exadg): High-Order Discontinuous Galerkin for the Exa-Scale
- [exapde/Exasim](https://github.com/exapde/Exasim): Generating Discontinuous Galerkin Codes For Extreme Scalable Simulations
- [flexi-framework/flexi](https://github.com/flexi-framework/flexi): Open Source High-Order Unstructured Discontinuous Galerkin Fluid Dynamics Solver
- [Trixi.jl](https://github.com/trixi-framework/Trixi.jl): Adaptive high-order numerical simulations of hyperbolic PDEs in Julia
- [ChiDG](https://github.com/nwukie/ChiDG): A Chimera-based, discontinuous Galerkin solver

#### Other techniques

- [jostbr/shallow-water](https://github.com/jostbr/shallow-water): Python model solving the shallow water equations (linear momentum, nonlinear continuity)
- [PavelDoGreat/WebGL-Fluid-Simulation](https://github.com/PavelDoGreat/WebGL-Fluid-Simulation): Play with fluids in your browser (works even on mobile)
- [PyFR/PyFR](https://github.com/PyFR/PyFR): Framework for solving advection-diffusion type problems on streaming architectures using the Flux Reconstruction approach of Huynh.
- [precice/precice](https://github.com/precice/precice): A coupling library for partitioned multi-physics simulations, including, but not restricted to fluid-structure interaction and conjugate heat transfer simulations.
- [cwrowley/ibpm](https://github.com/cwrowley/ibpm): Immersed Boundary Projection Method (IBPM)
- [barbagroup/PetIBM](https://github.com/barbagroup/PetIBM): PetIBM - toolbox and applications of the immersed-boundary method on distributed-memory architectures
- [vortexmethods/VM2D](https://github.com/vortexmethods/VM2D): VM2D: Open-Source Code for 2D Flow Simulation by Using Meshless Lagrangian Vortex Methods
- [markstock/vic2d](https://github.com/markstock/vic2d): Two-dimensional semi-Lagrangian vortex method for very low viscosity fluid simulation
- [Cantera/cantera](https://github.com/Cantera/cantera): Chemical kinetics, thermodynamics, and transport tool suite
- [NREL/EnergyPlus](https://github.com/NREL/EnergyPlus): EnergyPlus™ is a whole building energy simulation program that engineers, architects, and researchers use to model both energy consumption and water use in buildings.
- [uDALES/u-dales](https://github.com/uDALES/u-dales): uDALES: large-eddy-simulation software for urban flow, dispersion and microclimate modelling
- [taichi-dev/taichi](https://github.com/taichi-dev/taichi):  Parallel programming for everyone.
- [DARcorporation/xfoil-python](https://github.com/DARcorporation/xfoil-python): Stripped down version of XFOIL as compiled python module
- [mdolab/dafoam](https://github.com/mdolab/dafoam):  DAFoam: Discrete Adjoint with OpenFOAM for High-fidelity Gradient-based Design Optimization
- [peterdsharpe/AeroSandbox](https://github.com/peterdsharpe/AeroSandbox): Aircraft design optimization made fast through modern automatic differentiation. Plug-and-play analysis tools for aerodynamics, propulsion, structures, trajectory design, and much more.
- [team-ocean/veros](https://github.com/team-ocean/veros): The versatile ocean simulator, in pure Python, powered by JAX.
- [CliMA/Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl):
software for fast, friendly, flexible, data-driven, ocean-flavored fluid dynamics on CPUs and GPUs
- [gdeskos/DVMpp](https://github.com/gdeskos/DVMpp): 2D Discrete Vortex Method Code written in C++
- [enzo-project/enzo-dev](https://github.com/enzo-project/enzo-dev): The Enzo adaptive mesh-refinement simulation code
- [AMReX](https://amrex-codes.github.io/amrex/index.html): A software framework for massively parallel, block-structured adaptive mesh refinement (AMR) applications
- [llnl/samrai](https://github.com/llnl/samrai): Structured Adaptive Mesh Refinement Application Infrastructure - a scalable C++ framework for block-structured AMR application development
- [ForestClaw](http://www.forestclaw.org/ForestClaw/index.html): A parallel, adaptive library for logically Cartesian, mapped, multiblock domains
- [Hawkmoth flight simulation](https://biomech.web.unc.edu/hawkmoth-flight-simulation/): A blade-element model of a hawkmoth (Manduca sexta) implemented as an ODE in MATLAB
- [CFDEM@project](https://www.cfdem.com/): Dedicated to open source high performance scientific computing in fluid mechanics and particle science
- [adamslab-ub/optimizing-bio-inspired-riblets-openfoam](https://github.com/adamslab-ub/optimizing-bio-inspired-riblets-openfoam): This repository contains key software artifacts of the Open-FOAM based CFD framework that computes the drag coefficient of a NACA-0012 airfoil with bio-inspired surface riblets.
- [andrealani/COOLFluiD](https://github.com/andrealani/COOLFluiD): The object-oriented HPC platform for CFD, plasma and multi-physics simulations
- [Fluidity](https://fluidityproject.github.io/): An open-source computational fluid dynamics code with adaptive unstructured mesh capabilities
- [Basilisk](http://basilisk.fr/): Basilisk is the solution of partial differential equations on adaptive Cartesian meshes
- [Cantera/cantera](https://github.com/Cantera/cantera): Chemical kinetics, thermodynamics, and transport tool suite
- [code-mphi/ECOGEN](https://github.com/code-mphi/ECOGEN): A CFD open source code dedicated to multiphase compressible flows
- [Code_Aster](https://code-aster.org/spip.php?rubrique2): Structures and Thermomechanics Analysis for Studies and Research
- [GRAL](https://gral.tugraz.at/): The initial driver for the development of GRAL was the need for a model that could deal with the frequent low-wind-speed conditions in the inner-Alpine basins of Austria. Another important feature of GRAL is the ability to deal with the dispersion of pollutants emitted from road tunnel portals
- [fluiddyn/fluiddyn](https://foss.heptapod.net/fluiddyn/fluiddyn): FluidDyn project is an ecosystem of packages for research and teaching in fluid dynamics
- [fluiddyn/fluidsim](https://github.com/fluiddyn/fluidsim): Framework for studying fluid dynamics with numerical simulations using Python
- [gismo/gismo](https://github.com/gismo/gismo): A C++ library for isogeometric analysis
- [HiFiLES](https://hifiles.stanford.edu/): High Fidelity Large Eddy Simulation
- [HiSA](https://hisa.gitlab.io/): C++ based tool for computing compressible transonic and supersonic flow
- [FLOWUnsteady](https://flow.byu.edu/FLOWUnsteady/): A simulation engine of mixed-fidelity unsteady aerodynamics and aeroacoustics
- [REEF3D](https://reef3d.wordpress.com/): An open-source hydrodynamics framework
- [TurbGen](https://turbulence.utah.edu/code/): A collection of source code and tools for generating and postprocessing turbulence data
- [HPCWE](https://www.hpcwe-project.eu/codes/): High Performance Computing for Wind Energy
- [hyStrath](https://vincentcasseau.github.io/): Featuring hypersonic and rarefied gas dynamics
- [SHARPy](https://github.com/ImperialCollegeLondon/sharpy): Simulation of High Aspect Ratio aeroplanes and wind turbines in Python, a nonlinear aeroelastic code
- [Kratos](https://github.com/KratosMultiphysics/Kratos): Kratos Multiphysics is a framework for building parallel multi-disciplinary simulation software
- [KVSlab/turtleFSI](https://github.com/KVSlab/turtleFSI): Monolithic Fluid-Structure Interaction (FSI) solver
- [libmesh](https://libmesh.github.io/index.html): A framework for the numerical simulation of partial differential equations using arbitrary unstructured discretizations on serial and parallel platforms. A major goal of the library is to provide support for adaptive mesh refinement (AMR) computations in parallel while allowing a research scientist to focus on the physics they are modeling
- [MHYSA](https://github.com/LLNL/mhysa): Multispecies Hypersonic Flow Simulator
- [lucasschn/dstoolbox](https://github.com/lucasschn/dstoolbox): DSToolbox aims at creating a framework for the analysis of Dynamic Stall based on experimental data. The experimental data might be produced using an airfoil describing different type of motion, especially pitching and ramping up
- [Gerris](http://gfs.sourceforge.net/wiki/index.php/Main_Page): Gerris is a Free Software program for the solution of the partial differential equations describing fluid flow
- [matteobernardini/STREAmS](https://github.com/matteobernardini/STREAmS): Direct Numerical Simulations of compressible turbulent flows in Cartesian geometry solving the unsteady, fully compressible Navier-Stokes equations for a perfect gas
- [MBDyn](https://www.mbdyn.org/): General purpose Multibody Dynamics analysis software
- [mdolab/idwarp](https://github.com/mdolab/idwarp): IDWarp uses an inverse distance method to modify the location of mesh volume nodes given a perturbation of the surface nodes
- [MechSys](http://mechsys.nongnu.org/index.html): Library for the implementation of simulation tools in mechanics
- [MFC](https://mflowcode.github.io/): MFC is a fully-documented parallel simulation software for multi-component, multi-phase, and bubbly flows
- [MoDeNa-EUProject/MoDeNa](https://github.com/MoDeNa-EUProject/MoDeNa): About Software Framework for MOdelling of morphology DEvelopment of micro- and NAnostructures
- [momokhalil/pyHype](https://github.com/momokhalil/pyHype): Python framework for developing parallelized Computational Fluid Dynamics software to solve the hyperbolic 2D Euler equations on distributed, multi-block structured grids
- [moulin1024/WiRE-LES2](https://github.com/moulin1024/WiRE-LES2): Large-eddy simulation code written in CUDA Fortran for simulating atmospheric boundary layer flows
- [OPM Flow](https://opm-project.org/?page_id=19): OPM Flow is a fully-implicit, black-oil simulator capable of running industry-standard simulation models. The simulator is implemented using automatic differentiation to enable rapid development of new fluid models.
- [OuroPablo/Hydro3D](https://github.com/OuroPablo/Hydro3D): An open source Large Eddy Simulation code
- [PhysicsofFluids/AFiD](https://github.com/PhysicsofFluids/AFiD): A highly parallel application for Rayleigh-Benard and Taylor-Couette flows
- [smdogroup/tacs](https://github.com/smdogroup/tacs): parallel finite-element code for analysis and gradient-based design of advanced structures
- [UCNS3D](https://ucns3d.com/): Unstructured Grids, Compressible Solver, Navier–Stokes 3D
- [vectorfitting](https://www.sintef.no/projectweb/vectorfitting/): Vector Fitting is a robust numerical method for rational approximation in the frequency domain using poles and residues
- [uDALES/u-dales](https://github.com/uDALES/u-dales): Large-eddy-simulation software for urban flow, dispersion and microclimate modelling
- [vvflow/vvflow](https://github.com/vvflow/vvflow/): CFD software for performing flow simulations with Viscous Vortex Domains (VVD) method
- [weymouth/WaterLily.jl](https://github.com/weymouth/WaterLily.jl): Fast and simple fluid simulator in Julia
- [YsoSirius/windfarmGA](https://github.com/YsoSirius/windfarmGA): R Package to Optimize Windfarm Layouts
- [AMReX-Combustion/PeleC](https://github.com/AMReX-Combustion/PeleC):  an adaptive-mesh compressible hydrodynamics code for reacting flows
- [openCFS](https://opencfs.org/): a finite element-based multi-physics modelling and simulation tool
- [lukasumek/WRF_LES_diagnostics](https://github.com/lukasumek/WRF_LES_diagnostics/tree/master): modification of the WRF model source code (version 4.1) to derive turbulence and mean flow diagnostics
- [CFDEMproject/CFDEMcoupling-PUBLIC](https://github.com/CFDEMproject/CFDEMcoupling-PUBLIC): this code provides models and solvers to realize coupled CFD-DEM simulations using LIGGGHTS and OpenFOAM® technology
- [msolids/musen](https://github.com/msolids/musen): GPU-accelerated DEM simulation framework
- [quinoacomputing/quinoa](https://github.com/quinoacomputing/quinoa): Adaptive computational fluid dynamics
- [p-costa/CaNS](https://github.com/p-costa/CaNS): A code for fast, massively-parallel direct numerical simulations (DNS) of canonical flows
- [DNSLABIB](https://github.com/Aalto-CFD/DNSLABIB): DNSLABIB is a tool for Computational Fluid Dynamics (CFD), implementing the incompressible Navier-Stokes equations with a passive scalar transport equation and Lagrangian particle tracking in MATLAB both on a GPU and CPU
- [OpenHyperFLOW2D](https://github.com/sergeas67/OpenHyperFLOW2D): Parallel (C++/MPI/OpenMP/CUDA) research-educational CFD code for simulation 2D (flat/axisymmetrical) transient viscous compressible multicomponent sub/trans/supersonic reacting gas flow with RANS/URANS turbulence models
- [ECOGEN](https://github.com/code-mphi/ECOGEN): A CFD open source code dedicated to multiphase compressible flows
- [Mirheo](https://github.com/cselab/Mirheo): Mirheo is a GPU high-performance and high-throughput code aimed at simulation of flows at milli- and microscales
- [turtle](https://gitlab.mpcdf.mpg.de/TurTLE/turtle): TurTLE implements a number of standard functionality of Fourier-based pseudo-spectral numerical simulations, as well as the corresponding numerical particle-tracking functionality
- [pressel/pycles](https://github.com/pressel/pycles): A python based infrastructure for cloud large eddy simulation
- [mathLab/pi-BEM](https://github.com/mathLab/pi-BEM): Parallel Boundary Element Method solver
- [JULES](https://jules.jchmr.org/)
- [NGA2](https://github.com/jessecaps/NGA2): Object-oriented multi-mesh version of the classic reacting turbulent multiphase flow solver.

#### OpenFOAM Related
- [lesituu](https://bitbucket.org/lesituu/): Some libraries (wall model) for OpenFOAM
- [Caelus](http://www.caelus-cml.com/): Caelus is the next evolution in open-source computational continuum mechanics solutions. Caelus is a derivative of OpenFOAM® but has been restructured to create a stronger foundation on which to build your open-source CFD solutions
- [exaFOAM](https://exafoam.eu/): The ambitious exaFOAM project aims at overcoming the current limitations of Computational Fluid Dynamics (CFD) technology, especially in what concerns the exploitation of massively parallel HPC architectures.
- [mathLab/ITHACA-FV](https://github.com/mathLab/ITHACA-FV): Reduced order modelling techniques for OpenFOAM
- [BSL-EARSM](https://gitlab.com/ejp/BSL-EARSM/-/tree/master): The BSL-EARSM turbulence model of Menter etal (2012) (Menter, F. R., A. V. Garbaruk, and Y. Egorov. "Explicit algebraic Reynolds stress models for anisotropic wall-bounded flows." In Progress in Flight Physics, vol. 3, pp. 89-104. EDP Sciences, 2012.)
- [gmdpapers/hydrothermalfoam](https://gitlab.com/gmdpapers/hydrothermalfoam): Combination of hydrothermal and OpenFOAM — a three dimensional hydro-thermo-transport model designed to resolvefluid flow within submarine hydrothermal circulation systems
- [OpenQBMM](https://www.openqbmm.org/): A suite of solvers to simulate polydisperse multiphase flows using Quadrature-Based Moment Methods (QBMM) based on OpenFOAM
- [olaFlow](https://olaflow.github.io/): olaFlow CFD Suite is a free and open source project committed to bringing the latest advances for the simulation of wave dynamics to the OpenFOAM and FOAM-extend communities
- [OpenFOAM-ShihQuadraticKE](https://github.com/sagarsaroha18/OpenFOAM): code for ShihQuadraticKE turbulence model
- [Eddylicious ](https://github.com/timofeymukha/eddylicious): A python library for generating inflow boundary fields for scale-resolving simulations of turbulent flow like LES and DNS
- [Atizar/RapidCFD-dev](https://github.com/Atizar/RapidCFD-dev/tree/master): RapidCFD is an OpenFOAM fork running fully on CUDA platform
- [zhulianhua/dugksFoam](https://github.com/zhulianhua/dugksFoam): An OpenFOAM solver for Boltzmann model equation using discrete unified gas kinetic scheme
- [Kiiree/Adaptive-DES-model-2015](https://github.com/Kiiree/Adaptive-DES-model-2015): Source code for Yin, Zifei, K. R. Reddy, and Paul A. Durbin. "On the dynamic computation of the model constant in delayed detached eddy simulation." Physics of Fluids (1994-present) 27.2 (2015): 025105.
- [TonkomoLLC/chtMultiRegionReactingFoam](https://github.com/TonkomoLLC/chtMultiRegionReactingFoam): OpenFOAM transient solver for buoyant, turbulent fluid flow and solid heat conduction with conjugate heat transfer between solid and fluid regions, plus combustion with chemical reactions
- [ThomasKleiner/tpcMultiRegionFoam](https://github.com/ThomasKleiner/tpcMultiRegionFoam): CFD Solver for heat transfer simulations between solid and fluid regions with an implemented model for thermal phase change of pure substances in fluid regions. Algorithm is tested and validated for thin film applications with laminar flow.
- [Spationaute/SuperMarine](https://github.com/Spationaute/SuperMarine): Geometry Scipts for blockMesh in OpenFOAM
- [sedfoam/sedfoam](https://github.com/sedfoam/sedfoam): A three-dimensional two-phase flow solver, SedFoam, has been developed for sediment transport applications.
- [asimonder/projectionFoam](https://github.com/asimonder/projectionFoam): The pisoFoam solver to solve incompressible transient flows is iterative. Here it is replaced with an incremental projection scheme which is non-iterative therefore costs less than pisoFoam.
- [precice/openfoam-adapter](https://github.com/precice/openfoam-adapter): OpenFOAM-preCICE adapter
- [parallelwindfarms/cylinder](https://github.com/parallelwindfarms/cylinder): This package applies parallel-in-time integration, specifically Parareal, to OpenFOAM
- [mtezzele/MORtech-2019-Special-Issue](https://github.com/mtezzele/MORtech-2019-Special-Issue): All the material to reproduce the results of the work titled: "Enhancing CFD predictions in shape design problems by model and parameter space reduction".
- [krebeljk/openInjMoldSim](https://github.com/krebeljk/openInjMoldSim): Open source injection molding simulation. A solver for OpenFOAM.
- [joslorgom/shapeOptimizationFoam](https://github.com/joslorgom/shapeOptimizationFoam): Optimal Shape Design in External Flow with OpenFOAM.
- [joslorgom/heatAdjointFoam](https://github.com/joslorgom/heatAdjointFoam): Interior Control of the Heat Equation with the Steepest Descent Method in OpenFOAM.
- [hgopalan/RadianceToFoam](https://github.com/hgopalan/RadianceToFoam): Read the short-wave radiation information from Radiance and use it as a time-varying boundary condition in OpenFOAM.
- [isoAdvector/isoAdvector](https://github.com/isoAdvector/isoAdvector):  a geometric Volume-of-Fluid method for advection of a sharp interface between two incompressible fluids.
- [furstj/myTurbulenceModels](https://github.com/furstj/myTurbulenceModels): Additional turbulence models for OpenFOAM
- [Franjcf/hybridPorousInterFoam](https://github.com/Franjcf/hybridPorousInterFoam): OpenFOAM solver for performing single- and two-phase flow simulations on hybrid-scale porous media.
- [davidsblom/FOAM-FSI](https://github.com/davidsblom/FOAM-FSI): Fluid-Structure Interaction solvers for foam-extend
- [chenxijun1029/idugksFoam](https://github.com/chenxijun1029/idugksFoam): An OpenFOAM solver for incompressible isothermal fluid.
- [BlenderFOAM](https://nathanrooy.github.io/posts/2019-03-05/blenderfoam-aerodynamic-shape-optimization/): Open-source Fluid Based Shape Optimization
- [BjarkeEltardLarsen/stabRAS_OF50](https://github.com/BjarkeEltardLarsen/stabRAS_OF50): This repository contains formally stabilized versions of standard two-equation turbulence models for OpenFOAM-5.0
- [Rdfing/TurbulenceModel](https://github.com/Rdfing/TurbulenceModel): Some turbulence models used in OpenFOAM
- [classy_blocks](https://github.com/damogranlabs/classy_blocks): Python classes for easier creation of openFoam's blockMesh dictionaries.
- [blastfoam](https://github.com/synthetik-technologies/blastfoam): A CFD solver for multi-component compressible flow with application to high-explosive detonation, explosive safety and air blast.

##### Learn OpenFOAM

- [unicfdlab/TrainingTracks](https://github.com/unicfdlab/TrainingTracks): Materials for training tracks for continua media - OpenFOAM, vortex method, and other
- [phresher/OpenFOAM_Tutorials_Plus](https://github.com/phresher/OpenFOAM_Tutorials_Plus): OpenFOAM cases and notes
- [openfoamtutorials/openfoam_tutorials](https://github.com/openfoamtutorials/openfoam_tutorials): OpenFOAM tutorial cases with comprehensive instructions.
- [mrklein/foam-cases](https://github.com/mrklein/foam-cases): Cases created to illustrate certain OpenFOAM(R)/foam-extend functionality.
- [OS_CFD](http://www.tfd.chalmers.se/~hani/kurser/OS_CFD/): PhD course CFD with OpenSource Software

### GPU

- [Fiesta](https://zenodo.org/record/4554565): FIESTA (Fast InterfacES and Transition in the Atmosphere) is a GPU accelerated, exascale ready computational fluid dynamics code.  FIESTA is designed to solve compressible flow problems with gas interfaces at a range of physical scales.
- [HTR-solver](https://github.com/stanfordhpccenter/HTR-solver): Hypersonic Task-based Research (HTR) solver for the Navier-Stokes equations at hypersonic Mach numbers including finite-rate chemistry for dissociating air and multicomponent transport.
- [maxcuda/CaNS](https://github.com/maxcuda/CaNS): A code for fast, massively-parallel direct numerical simulations (DNS) of canonical flows
- [Nek5000/nekRS](https://github.com/Nek5000/nekRS): Navier Stokes solver based on the spectral element method targeting classical processors and hardware accelerators like GPUs
- [AFiD_GPU_opensource](https://github.com/PhysicsofFluids/AFiD_GPU_opensource): Version of the AFiD code ported to GPU accelerators


### Rotor
- [dust](https://gitlab.com/dust_group/dust): a flexible solution for aerodynamic simulation of complex configurations

### Machine learning / Deep Learning

- [lululxvi/deepxde](https://github.com/lululxvi/deepxde): Deep learning library for solving differential equations and more.
- [SciML/NeuralPDE.jl](https://github.com/SciML/NeuralPDE.jl): Physics-Informed Neural Networks (PINN) and Deep BSDE Solvers of Differential Equations for Scientific Machine Learning (SciML) accelerated simulation.
- [google/neural-tangents](https://github.com/google/neural-tangents):  Fast and Easy Infinite Neural Networks in Python.
- [cornellius-gp/gpytorch](https://github.com/cornellius-gp/gpytorch): A highly efficient and modular implementation of Gaussian Processes in PyTorch.
- [google/jax-cfd](https://github.com/google/jax-cfd): Computational Fluid Dynamics in JAX.
- [isl-org/DeepLagrangianFluids](https://github.com/isl-org/DeepLagrangianFluids): Lagrangian Fluid Simulation with Continuous Convolutions.
- [tum-pbs/PhiFlow](https://github.com/tum-pbs/PhiFlow): A differentiable PDE solving framework for machine learning.
- [maxjiang93/space_time_pde](https://github.com/maxjiang93/space_time_pde): MeshfreeFlowNet: Physical Constrained Space Time Super-Resolution
- [JAXFLUIDS](https://github.com/tumaer/JAXFLUIDS): Differentiable Fluid Dynamics Package.

#### Super-resolution

- [Jianxun-Wang/PICNNSR](https://github.com/Jianxun-Wang/PICNNSR): About Super-resolution and denoising of fluid flow using physics-informed convolutional neural networks without high-resolution labels -- parametric forward SR and boundary inference

#### Uncertainty Quantification

- [DAKOTA](https://dakota.sandia.gov/): The Dakota project delivers both state-of-the-art research and robust, usable software for optimization and UQ
- [multiUQ](https://bitbucket.org/markowkes/multiuq/src/master/): Multiphase flow solver with uncertainty quantification
- [KTH-Nek5000/UQit](https://github.com/KTH-Nek5000/UQit): A Python Package for Uncertainty Quantification (UQ) in Computational Fluid Dynamics
- [SURGroup/UQpy](https://github.com/SURGroup/UQpy): UQpy (Uncertainty Quantification with python) is a general purpose Python toolbox for modeling uncertainty in physical and mathematical systems.

#### Optimization

- [OpenMDAO](http://openmdao.org): OpenMDAO is an open-source high-performance computing platform for systems analysis and multidisciplinary optimization, written in Python

#### Platforms

- [kokkos](https://github.com/kokkos/kokkos): Kokkos C++ Performance Portability Programming EcoSystem: The Programming Model - Parallel Execution and Memory Abstraction

## 4. Tools

- [mathLab/PyDMD](https://github.com/mathLab/PyDMD): Python Dynamic Mode Decomposition
- [dynamicslab/pysindy](https://github.com/dynamicslab/pysindy): A sparse regression package with several implementations for the Sparse Identification of Nonlinear Dynamical systems.
- [ritchieng/eigenvectors-from-eigenvalues](https://github.com/ritchieng/eigenvectors-from-eigenvalues): This repository implements a calculation of eigenvectors from eigenvectors elegantly through PyTorch.
- [mengaldo/PySPOD](https://github.com/mengaldo/PySPOD): A Python package for spectral proper orthogonal decomposition (SPOD).
- [belson17/modred](https://github.com/belson17/modred): An easy-to-use and parallelized library for finding modal decompositions and reduced-order models.
- [Astroua/TurbuStat](https://github.com/Astroua/TurbuStat): Statistics of Turbulence Python Package.
- [haller-group/LCStool](https://github.com/haller-group/LCStool): LCStool: LCS Tool is a computational engine for analyzing fluid flows by extracting their most influential material surfaces, Lagrangian Coherent Structures.
- [mauritssilvis/lesTools](https://github.com/mauritssilvis/lesTools): A toolbox for the construction and assessment of subgrid-scale models for large-eddy simulations
- [AtzoriMarco/InSituPackage](https://github.com/AtzoriMarco/InSituPackage): In-Situ analysis for Nek5000 (v17) and Catalyst
- [GIBBON](https://www.gibboncode.org/): IBBON (The Geometry and Image-Based Bioengineering add-On) is an open-source MATLAB toolbox by Kevin M. Moerman and includes an array of image and geometry visualization and processing tools and is interfaced with free open source software such as TetGen, for robust tetrahedral meshing, and FEBio for finite element analysis
- [fluiddyn/fluidfoam](https://github.com/fluiddyn/fluidfoam): OpenFoam postprocessing python tool
- [inendi/inspector](https://gitlab.com/inendi/inspector): INENDI Inspector allows any user to perform interactive explorations over very large amounts of data.
- [redbKIT](http://redbkit.github.io/redbKIT/): A MATLAB library for reduced-order modeling of PDEs.

## 6. Datasets

- [shengzesnail/PIV_dataset](https://github.com/shengzesnail/PIV_dataset): PIV dataset
- [Johns Hopkins Turbulence Databases](http://turbulence.pha.jhu.edu/)
- [idies/pyJHTDB](https://github.com/idies/pyJHTDB): Python wrapper for the Johns Hopkins turbulence database library.
- [Turbulence Database Madrid](https://torroja.dmt.upm.es/turbdata/): Open turbulence dataset provided by Fluid Dynamics Group UPM.
- [Kawamura Lab - DNS Database of Wall Turbulence and Heat Transfer](https://www.rs.tus.ac.jp/t2lab/db/index.html)
- [Turbulence Database UTexas](https://turbulence.oden.utexas.edu/)
- [Boundary Layer DNS/LES Data KTH](https://www.mech.kth.se/~pschlatt/DATA/)
- [DNS database of turbulent flows](http://newton.dma.uniroma1.it/)
- [Turbulence modeling gateway](http://turbgate.engin.umich.edu/)

## 8. Related links

- [alexlib/awesome_piv](https://github.com/alexlib/awesome_piv); A curated list of repositories related to PIV (particle image velocimetry).
- [nschloe/awesome-scientific-computing](https://github.com/nschloe/awesome-scientific-computing): Curated list of awesome software for numerical analysis and scientific computing.
- [qd-cae/awesome-CAE](https://github.com/qd-cae/awesome-CAE): A curated list of awesome CAE frameworks, libraries and software.
- [awesomedata/awesome-public-datasets](https://github.com/awesomedata/awesome-public-datasets):  A topic-centric list of HQ open datasets.
- [Open source software for CAE](https://www.cfdsupport.com/cae-open-source-software.html): Open source software for CAE.
- [turbulence-codes](https://github.com/sthavishtha/turbulence-codes): Curated list of some open-source codes for turbulent flow simulations, including turbulent multiphase, turbulent reacting flows, turbulent convection and turbulent atmospheric physics.
- [SPH CODES](https://www.spheric-sph.org/sph-projects-and-codes): SPH code & software
- [CFD Machine Learning Useful Links](https://simulationwork.com/2021/12/23/cfd-machine-learning-useful-links/): A curated list of machine learning papers, codes, libraries, and databases applied to fluid mechanics.
- [Center for Engineering and Scientific Computation](http://www.cesc.zju.edu.cn/learningcenter/fluidlist_e.htm)
- [CFD Machine Learning Useful Links](https://simulationwork.com/2021/12/23/cfd-machine-learning-useful-links/)
- [Machine learning for turbulence modeling](https://www.monolithai.com/post/machine-learning-for-turbulence-modeling)
- [Open-source Computational Physics](https://holger-marschall.owncube.com/blog/2022/07/04/open-source-computational-physics-new-kids-on-the-block-and-the-usual-suspects/)
- [Open-source code shared by Prof. Fu Lin](http://linfu.people.ust.hk/open-source-codes.html)
- [Open-source code shared by CTFLab](https://www.multiphasecfd.com/codes)
- [lento234/awesome-fluid-dynamics](https://github.com/lento234/awesome-fluid-dynamics): A curated list of repositories related to fluid dynamics.
