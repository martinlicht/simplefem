SHELL = /bin/sh

default: all

include ../common.recipe.mk 

Objects.basic = 
Objects.combinatorics = ../combinatorics/*.o
Objects.operators = ../operators/*.o
Objects.dense = ../dense/*.o
Objects.sparse = ../sparse/*.o
Objects.solver = ../solver/*.o
Objects.mesh = ../mesh/*.o
Objects.vtk = ../vtk/*.o
Objects.matrixmarket = ../matrixmarket/*.o
Objects.fem =  

Links.basic         = $(Objects.basic) 
Links.combinatorics = $(Objects.combinatorics) 
Links.operators     = $(Objects.operators) $(Objects.combinatorics)
Links.dense         = $(Objects.dense) $(Objects.sparse) $(Objects.operators) $(Objects.combinatorics)
Links.sparse        = $(Objects.sparse) $(Objects.operators) $(Objects.combinatorics) $(Objects.dense)
Links.solver        = $(Objects.solver) $(Objects.operators) $(Objects.sparse) $(Objects.dense) $(Objects.combinatorics)
Links.mesh          = $(Objects.mesh) $(Objects.combinatorics) $(Objects.operators) $(Objects.dense) $(Objects.sparse) 
Links.vtk           = $(Objects.vtk) $(Objects.mesh) $(Objects.dense) $(Objects.sparse) $(Objects.combinatorics) $(Objects.operators)
Links.matrixmarket  = $(Objects.matrixmarket) 
Links.fem           = $(Objects.fem) $(Objects.vtk) $(Objects.mesh) $(Objects.solver) $(Objects.dense) $(Objects.sparse) $(Objects.combinatorics) $(Objects.operators)


./*/*.out: makefile ../common.recipe.mk


depdir := .deps





basic.sources := $(wildcard basic/*.cpp)

basic.outs := $(patsubst %.cpp,%.out,$(basic.sources))

basic.depdir := basic/$(depdir)

basic.dependencies := $(patsubst basic/%.cpp,basic/.deps/%.d,$(basic.sources))

$(basic.depdir): ; @mkdir -p $@

$(basic.outs): basic/%.out: basic/%.cpp | $(basic.depdir)
	g++ -MM basic/$*.cpp -MT $@ -MF basic/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.basic) -o $@ 

-include $(basic.dependencies)

basic.tests: $(basic.outs)





combinatorics.sources := $(wildcard combinatorics/*.cpp)

combinatorics.outs := $(patsubst %.cpp,%.out,$(combinatorics.sources))

combinatorics.depdir := combinatorics/$(depdir)

combinatorics.dependencies := $(patsubst combinatorics/%.cpp,combinatorics/.deps/%.d,$(combinatorics.sources))

$(combinatorics.depdir): ; @mkdir -p $@

$(combinatorics.outs): combinatorics/%.out: combinatorics/%.cpp | $(combinatorics.depdir)
	@g++ -MM combinatorics/$*.cpp -MT $@ -MF combinatorics/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.combinatorics) -o $@ 

-include $(combinatorics.dependencies)

combinatorics.tests: $(combinatorics.outs)





    



operators.sources := $(wildcard operators/*.cpp)

operators.outs := $(patsubst %.cpp,%.out,$(operators.sources))

operators.depdir := operators/$(depdir)

operators.dependencies := $(patsubst operators/%.cpp,operators/.deps/%.d,$(operators.sources))

$(operators.depdir): ; @mkdir -p $@

$(operators.outs): operators/%.out: operators/%.cpp | $(operators.depdir)
	@g++ -MM operators/$*.cpp -MT $@ -MF operators/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.operators) -o $@ 

-include $(operators.dependencies)

operators.tests: $(operators.outs)



    
    
    
    
dense.sources := $(wildcard dense/*.cpp)

dense.outs := $(patsubst %.cpp,%.out,$(dense.sources))

dense.depdir := dense/$(depdir)

dense.dependencies := $(patsubst dense/%.cpp,dense/.deps/%.d,$(dense.sources))

$(dense.depdir): ; @mkdir -p $@

$(dense.outs): dense/%.out: dense/%.cpp | $(dense.depdir)
	@g++ -MM dense/$*.cpp -MT $@ -MF dense/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.dense) -o $@ 

-include $(dense.dependencies)

dense.tests: $(dense.outs)








sparse.sources := $(wildcard sparse/*.cpp)

sparse.outs := $(patsubst %.cpp,%.out,$(sparse.sources))

sparse.depdir := sparse/$(depdir)

sparse.dependencies := $(patsubst sparse/%.cpp,sparse/.deps/%.d,$(sparse.sources))

$(sparse.depdir): ; @mkdir -p $@

$(sparse.outs): sparse/%.out: sparse/%.cpp | $(sparse.depdir)
	@g++ -MM sparse/$*.cpp -MT $@ -MF sparse/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.sparse) -o $@ 

-include $(sparse.dependencies)

sparse.tests: $(sparse.outs)






solver.sources := $(wildcard solver/*.cpp)

solver.outs := $(patsubst %.cpp,%.out,$(solver.sources))

solver.depdir := solver/$(depdir)

solver.dependencies := $(patsubst solver/%.cpp,solver/.deps/%.d,$(solver.sources))

$(solver.depdir): ; @mkdir -p $@

$(solver.outs): solver/%.out: solver/%.cpp | $(solver.depdir)
	@g++ -MM solver/$*.cpp -MT $@ -MF solver/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.solver) -o $@ 

-include $(solver.dependencies)

solver.tests: $(solver.outs)






mesh.sources := $(wildcard mesh/*.cpp)

mesh.outs := $(patsubst %.cpp,%.out,$(mesh.sources))

mesh.depdir := mesh/$(depdir)

mesh.dependencies := $(patsubst mesh/%.cpp,mesh/.deps/%.d,$(mesh.sources))

$(mesh.depdir): ; @mkdir -p $@

$(mesh.outs): mesh/%.out: mesh/%.cpp | $(mesh.depdir)
	@g++ -MM mesh/$*.cpp -MT $@ -MF mesh/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.mesh) -o $@ 

-include $(mesh.dependencies)

mesh.tests: $(mesh.outs)






vtk.sources := $(wildcard vtk/*.cpp)

vtk.outs := $(patsubst %.cpp,%.out,$(vtk.sources))

vtk.depdir := vtk/$(depdir)

vtk.dependencies := $(patsubst vtk/%.cpp,vtk/.deps/%.d,$(vtk.sources))

$(vtk.depdir): ; @mkdir -p $@

$(vtk.outs): vtk/%.out: vtk/%.cpp | $(vtk.depdir)
	@g++ -MM vtk/$*.cpp -MT $@ -MF vtk/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.vtk) -o $@ 

-include $(vtk.dependencies)

vtk.tests: $(vtk.outs)






matrixmarket.sources := $(wildcard matrixmarket/*.cpp)

matrixmarket.outs := $(patsubst %.cpp,%.out,$(matrixmarket.sources))

matrixmarket.depdir := matrixmarket/$(depdir)

matrixmarket.dependencies := $(patsubst matrixmarket/%.cpp,matrixmarket/.deps/%.d,$(matrixmarket.sources))

$(matrixmarket.depdir): ; @mkdir -p $@

$(matrixmarket.outs): matrixmarket/%.out: matrixmarket/%.cpp | $(matrixmarket.depdir)
	@g++ -MM matrixmarket/$*.cpp -MT $@ -MF matrixmarket/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.matrixmarket) -o $@ 

-include $(matrixmarket.dependencies)

matrixmarket.tests: $(matrixmarket.outs)





fem.sources := $(wildcard fem/*.cpp)

fem.outs := $(patsubst %.cpp,%.out,$(fem.sources))

fem.depdir := fem/$(depdir)

fem.dependencies := $(patsubst fem/%.cpp,fem/.deps/%.d,$(fem.sources))

$(fem.depdir): ; @mkdir -p $@

$(fem.outs): fem/%.out: fem/%.cpp | $(fem.depdir)
	@g++ -MM fem/$*.cpp -MT $@ -MF fem/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.fem) -o $@ 

-include $(fem.dependencies)

fem.tests: $(fem.outs)






# test.compiler.out: test.compiler.cpp 
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test.compiler.cpp  -o test.compiler.out
# 
# test.container.out: test.container.cpp 
# 	$(CXX) $(CXXFLAGS) $(CPPFLAGS) test.container.cpp  -o test.container.out



all: basic.tests 
all: combinatorics.tests
all: operators.tests
all: dense.tests
all: sparse.tests
all: solver.tests
all: mesh.tests
all: vtk.tests
all: matrixmarket.tests
all: fem.tests
# all: test.compiler.out 


include makefile.tests


clean: 
	-rm -f *.o *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 
	-rm -f */*.o */*.gch
	-rm -f */*.exe */*.exe.stackdump
	-rm -f */*.out */*.out.stackdump 


