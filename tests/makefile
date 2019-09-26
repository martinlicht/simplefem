SHELL = /bin/sh

default: all

include ../common.recipe.mk 

affix.combinatorics = combinatorics
affix.operators := operators combinatorics
affix.dense := dense sparse operators combinatorics
affix.sparse := sparse dense operators combinatorics
affix.solver := solver sparse dense operators combinatorics
affix.mesh := mesh dense sparse operators combinatorics
affix.vtk := vtk $(affix.mesh)
affix.matrixmarket = matrixmarket
affix.fem := vtk mesh solver dense sparse operators combinatorics

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

combinatorics.include := $(patsubst %,-L../%,$(affix.combinatorics))
linkerprefix=-Wl,
combinatorics.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.combinatorics)) 
combinatorics.rpath := $(patsubst %,$(linkerprefix)%,$(combinatorics.rpath_t)) 
combinatorics.link := $(patsubst %,-l%,$(affix.combinatorics))

$(combinatorics.depdir): ; @mkdir -p $@

$(combinatorics.outs): combinatorics/%.out: combinatorics/%.cpp | $(combinatorics.depdir)
	# @echo $(combinatorics.include)
	# @echo $(combinatorics.rpath)
	# @echo $(combinatorics.link)
	@g++ -MM combinatorics/$*.cpp -MT $@ -MF combinatorics/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.combinatorics) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(combinatorics.include) $(combinatorics.rpath) $(combinatorics.link) -o $@

-include $(combinatorics.dependencies)

combinatorics.tests: $(combinatorics.outs)





    



operators.sources := $(wildcard operators/*.cpp)

operators.outs := $(patsubst %.cpp,%.out,$(operators.sources))

operators.depdir := operators/$(depdir)

operators.dependencies := $(patsubst operators/%.cpp,operators/.deps/%.d,$(operators.sources))

operators.include := $(patsubst %,-L../%,$(affix.operators))
linkerprefix=-Wl,
operators.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.operators)) 
operators.rpath := $(patsubst %,$(linkerprefix)%,$(operators.rpath_t)) 
operators.link := $(patsubst %,-l%,$(affix.operators))

$(operators.depdir): ; @mkdir -p $@

$(operators.outs): operators/%.out: operators/%.cpp | $(operators.depdir)
	# @echo $(operators.include)
	# @echo $(operators.rpath)
	# @echo $(operators.link)
	@g++ -MM operators/$*.cpp -MT $@ -MF operators/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.operators) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(operators.include) $(operators.rpath) $(operators.link) -o $@

-include $(operators.dependencies)

operators.tests: $(operators.outs)



    
    
    
    
dense.sources := $(wildcard dense/*.cpp)

dense.outs := $(patsubst %.cpp,%.out,$(dense.sources))

dense.depdir := dense/$(depdir)

dense.dependencies := $(patsubst dense/%.cpp,dense/.deps/%.d,$(dense.sources))

dense.include := $(patsubst %,-L../%,$(affix.dense))
linkerprefix=-Wl,
dense.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.dense)) 
dense.rpath := $(patsubst %,$(linkerprefix)%,$(dense.rpath_t)) 
dense.link := $(patsubst %,-l%,$(affix.dense))

$(dense.depdir): ; @mkdir -p $@

$(dense.outs): dense/%.out: dense/%.cpp | $(dense.depdir)
	# @echo $(dense.include)
	# @echo $(dense.rpath)
	# @echo $(dense.link)
	@g++ -MM dense/$*.cpp -MT $@ -MF dense/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.dense) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(dense.include) $(dense.rpath) $(dense.link) -o $@

-include $(dense.dependencies)

dense.tests: $(dense.outs)








sparse.sources := $(wildcard sparse/*.cpp)

sparse.outs := $(patsubst %.cpp,%.out,$(sparse.sources))

sparse.depdir := sparse/$(depdir)

sparse.dependencies := $(patsubst sparse/%.cpp,sparse/.deps/%.d,$(sparse.sources))

sparse.include := $(patsubst %,-L../%,$(affix.sparse))
linkerprefix=-Wl,
sparse.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.sparse)) 
sparse.rpath := $(patsubst %,$(linkerprefix)%,$(sparse.rpath_t)) 
sparse.link := $(patsubst %,-l%,$(affix.sparse))

$(sparse.depdir): ; @mkdir -p $@

$(sparse.outs): sparse/%.out: sparse/%.cpp | $(sparse.depdir)
	# @echo $(sparse.include)
	# @echo $(sparse.rpath)
	# @echo $(sparse.link)
	@g++ -MM sparse/$*.cpp -MT $@ -MF sparse/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.sparse) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(sparse.include) $(sparse.rpath) $(sparse.link) -o $@

-include $(sparse.dependencies)

sparse.tests: $(sparse.outs)






solver.sources := $(wildcard solver/*.cpp)

solver.outs := $(patsubst %.cpp,%.out,$(solver.sources))

solver.depdir := solver/$(depdir)

solver.dependencies := $(patsubst solver/%.cpp,solver/.deps/%.d,$(solver.sources))

solver.include := $(patsubst %,-L../%,$(affix.solver))
linkerprefix=-Wl,
solver.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.solver)) 
solver.rpath := $(patsubst %,$(linkerprefix)%,$(solver.rpath_t)) 
solver.link := $(patsubst %,-l%,$(affix.solver))

$(solver.depdir): ; @mkdir -p $@

$(solver.outs): solver/%.out: solver/%.cpp | $(solver.depdir)
	# @echo $(solver.include)
	# @echo $(solver.rpath)
	# @echo $(solver.link)
	@g++ -MM solver/$*.cpp -MT $@ -MF solver/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.solver) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(solver.include) $(solver.rpath) $(solver.link) -o $@

-include $(solver.dependencies)

solver.tests: $(solver.outs)






mesh.sources := $(wildcard mesh/*.cpp)

mesh.outs := $(patsubst %.cpp,%.out,$(mesh.sources))

mesh.depdir := mesh/$(depdir)

mesh.dependencies := $(patsubst mesh/%.cpp,mesh/.deps/%.d,$(mesh.sources))

mesh.include := $(patsubst %,-L../%,$(affix.mesh))
linkerprefix=-Wl,
mesh.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.mesh)) 
mesh.rpath := $(patsubst %,$(linkerprefix)%,$(mesh.rpath_t)) 
mesh.link := $(patsubst %,-l%,$(affix.mesh))

$(mesh.depdir): ; @mkdir -p $@

$(mesh.outs): mesh/%.out: mesh/%.cpp | $(mesh.depdir)
	# @echo $(mesh.include)
	# @echo $(mesh.rpath)
	# @echo $(mesh.link)
	@g++ -MM mesh/$*.cpp -MT $@ -MF mesh/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.mesh) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(mesh.include) $(mesh.rpath) $(mesh.link) -o $@

-include $(mesh.dependencies)

mesh.tests: $(mesh.outs)






vtk.sources := $(wildcard vtk/*.cpp)

vtk.outs := $(patsubst %.cpp,%.out,$(vtk.sources))

vtk.depdir := vtk/$(depdir)

vtk.dependencies := $(patsubst vtk/%.cpp,vtk/.deps/%.d,$(vtk.sources))

vtk.include := $(patsubst %,-L../%,$(affix.vtk))
linkerprefix=-Wl,
vtk.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.vtk)) 
vtk.rpath := $(patsubst %,$(linkerprefix)%,$(vtk.rpath_t)) 
vtk.link := $(patsubst %,-l%,$(affix.vtk))

$(vtk.depdir): ; @mkdir -p $@

$(vtk.outs): vtk/%.out: vtk/%.cpp | $(vtk.depdir)
	# @echo $(vtk.include)
	# @echo $(vtk.rpath)
	# @echo $(vtk.link)
	@g++ -MM vtk/$*.cpp -MT $@ -MF vtk/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.vtk) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(vtk.include) $(vtk.rpath) $(vtk.link) -o $@

-include $(vtk.dependencies)

vtk.tests: $(vtk.outs)






matrixmarket.sources := $(wildcard matrixmarket/*.cpp)

matrixmarket.outs := $(patsubst %.cpp,%.out,$(matrixmarket.sources))

matrixmarket.depdir := matrixmarket/$(depdir)

matrixmarket.dependencies := $(patsubst matrixmarket/%.cpp,matrixmarket/.deps/%.d,$(matrixmarket.sources))

matrixmarket.include := $(patsubst %,-L../%,$(affix.matrixmarket))
linkerprefix=-Wl,
matrixmarket.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.matrixmarket)) 
matrixmarket.rpath := $(patsubst %,$(linkerprefix)%,$(matrixmarket.rpath_t)) 
matrixmarket.link := $(patsubst %,-l%,$(affix.matrixmarket))

$(matrixmarket.depdir): ; @mkdir -p $@

$(matrixmarket.outs): matrixmarket/%.out: matrixmarket/%.cpp | $(matrixmarket.depdir)
	# @echo $(matrixmarket.include)
	# @echo $(matrixmarket.rpath)
	# @echo $(matrixmarket.link)
	@g++ -MM matrixmarket/$*.cpp -MT $@ -MF matrixmarket/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.matrixmarket) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(matrixmarket.include) $(matrixmarket.rpath) $(matrixmarket.link) -o $@

-include $(matrixmarket.dependencies)

matrixmarket.tests: $(matrixmarket.outs)





fem.sources := $(wildcard fem/*.cpp)

fem.outs := $(patsubst %.cpp,%.out,$(fem.sources))

fem.depdir := fem/$(depdir)

fem.dependencies := $(patsubst fem/%.cpp,fem/.deps/%.d,$(fem.sources))

fem.include := $(patsubst %,-L../%,$(affix.fem))
linkerprefix:=-Wl,
fem.rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.fem)) 
fem.rpath := $(patsubst %,$(linkerprefix)%,$(fem.rpath_t)) 
fem.link := $(patsubst %,-l%,$(affix.fem))

$(fem.depdir): ; @mkdir -p $@

$(fem.outs): fem/%.out: fem/%.cpp | $(fem.depdir)
	# @echo $(fem.include)
	@echo $(fem.rpath)
	# @echo $(fem.link)
	@g++ -MM fem/$*.cpp -MT $@ -MF fem/.deps/$*.d
	# $(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(Links.fem) -o $@ 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(fem.include) $(fem.rpath) $(fem.link) -o $@

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


include makefile.runtests


clean: 
	-rm -f *.o *.gch
	-rm -f *.exe *.exe.stackdump
	-rm -f *.out *.out.stackdump 
	-rm -f */.deps/*.d
	-rm -f */*.o */*.gch
	-rm -f */*.exe */*.exe.stackdump
	-rm -f */*.out */*.out.stackdump 


