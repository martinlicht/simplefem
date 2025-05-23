# For each directory of tests, we define a list of modules necessary for linking
# the executables. These also enter the dependencies of the makefile targets.
# Each new category of tests must be added here together with all necessary modules.

affix.external      := 
affix.base          := base
affix.utility       := utility base external
affix.combinatorics := combinatorics base external
affix.operators     := operators combinatorics utility base external
affix.dense         := dense sparse operators combinatorics utility base external
affix.sparse        := sparse dense operators combinatorics utility base external
affix.solver        := solver sparse dense operators combinatorics utility base external
affix.mesh          := mesh dense sparse operators combinatorics utility base external
affix.vtk           := vtk $(affix.mesh)
#affix.matrixmarket  := matrixmarket
affix.fem           := fem mesh solver dense sparse operators combinatorics utility base external
affix.image         := fem mesh solver dense sparse operators combinatorics utility base external
affix.solverfem     := fem mesh solver dense sparse operators combinatorics utility base external
affix.sullivan2D    := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.sullivan3D    := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.whitney2D     := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.whitney3D     := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.eigenvalue    := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.nullspace     := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.afem          := fem vtk mesh solver dense sparse operators combinatorics utility base external
affix.benchmark     := fem vtk mesh solver dense sparse operators combinatorics utility base external
