# # Example usage:
# # 
# # 
# # include ../../common.compile.mk 
# # 
# # include ../../common.upkeep.mk
# # 
# # depdir := .deps
# # 
# # context=dense
# # contextdir=.
# # affix.$(context)=dense sparse operators combinatorics
# # projectdir=../../
# # pathvar=$(shell pwd)/../../

ifeq ($(OS),Windows_NT)
ending := exe
else
ending := out
endif

$(context).sources := $(sort $(wildcard $(contextdir)/*.cpp))

$(context).outs    := $(patsubst %.cpp,%.$(ending),$($(context).sources))

$(context).depdir  := $(contextdir)/$(depdir)

$(context).dependencies := $(patsubst $(contextdir)/%.cpp,$(contextdir)/.deps/%.d,$($(context).sources))

$(context).include := $(patsubst %,-L$(projectdir)/%,$(affix.$(context)))
linkerprefix       :=-Wl,
$(context).rpath_t := $(patsubst %,-rpath=$(pathvar)/%,$(affix.$(context))) 
$(context).rpath   := $(patsubst %,$(linkerprefix)%,$($(context).rpath_t)) 
$(context).lib     := $(patsubst %,-l%,$(affix.$(context)))
$(context).alib    := $(patsubst %,-l:lib%.a,$(affix.$(context)))
$(context).solib   := $(patsubst %,-l:lib%.so,$(affix.$(context)))

ifeq ($(LINKINGTYPE),static)
$(context).mylib := $($(context).alib)
else ifeq ($(LINKINGTYPE),dynamic)
$(context).mylib := $($(context).solib)
else ifeq ($(LINKINGTYPE),unspecified)
$(context).mylib := $($(context).lib)
else
$(error No linking mode recognized: $(LINKINGTYPE)) 
endif


.PHONY: $($(context).depdir)
$($(context).depdir):
	@-mkdir -p $@






$($(context).outs): mycontext    := $(context)
$($(context).outs): mycontextdir := $(contextdir)
$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/%.cpp | $($(context).depdir)
#	@ echo link type:   $(LINKINGTYPE)
#	@ echo context:     $(mycontext)
#	@ echo context dir: $(mycontextdir)
#	@ echo target:      $@
#	@ echo target:      $<
#	@ echo target:      $^
#	@ echo contextdir:  $(mycontextdir)
#	@ echo depdir:      $($(mycontext).depdir)
#	@ echo include:     $($(mycontext).include)
#	@ echo rpath:       $($(mycontext).rpath)
#	@ echo lib:         $($(mycontext).lib)
	@echo Compiling $@ ...
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -std=c++2a -MM $(mycontextdir)/$*.cpp -MT $@ -MF $($(mycontext).depdir)/$*.d
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include) $($(mycontext).rpath) $($(mycontext).mylib) -o $@ $(LDLIBS)
else
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include)                       $($(mycontext).mylib) -o $@ $(LDLIBS)
endif
#	cp $@ $(patsubst %.out,%.exe,$@)

-include $($(context).dependencies)

$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/*.mk $(projectdir)/tests/*.mk

$(context).tests: $($(context).outs)











$(context).sources     := $(sort $(wildcard $(contextdir)/*cpp))
$(context).executables := $(patsubst %.cpp,%.out,$($(context).sources))
$(context).runs        := $(patsubst %.cpp,%.run,$($(context).sources))
$(context).silent_runs := $(patsubst %.cpp,%.silent_run,$($(context).sources))


run: $(context).run

$(context).run: $($(context).runs)

$($(context).runs): %.run : %.out
	./$< 

silent_run: $(context).silent_run

$(context).silent_run: $($(context).silent_runs)

$($(context).silent_runs): %.silent_run : %.out
	./$< > /dev/null 

# 2> /dev/null

PHONY: $($(context).runs) $($(context).silent_runs)
