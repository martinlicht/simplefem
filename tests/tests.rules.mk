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
# cp $@ $(patsubst %.out,%.exe,$@)

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



#####################################
# How to compile the executables 

.PHONY: $($(context).depdir)
$($(context).depdir):
	@-mkdir -p $@

DEPFLAGS = -MT $@ -MF $($(mycontext).depdir)/$*.d -MP -MMD

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
#	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -std=c++2a -MT $@ -MF $($(mycontext).depdir)/$*.d -MM $(mycontextdir)/$*.cpp
ifeq ($(LINKINGTYPE),dynamic)
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include) $($(mycontext).rpath) $($(mycontext).mylib) -o $@ $(LDLIBS) $(DEPFLAGS)
else
	@$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include)                       $($(mycontext).mylib) -o $@ $(LDLIBS) $(DEPFLAGS)
endif
	
$($(context).dependencies):

-include $($(context).dependencies)

$($(context).outs): $(contextdir)/%.$(ending): $(contextdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/makefile 
$($(context).outs): $(contextdir)/%.$(ending): $(projectdir)/*.mk $(projectdir)/tests/*.mk

$(context).tests: $($(context).outs)





###################################################
# How to automatically run the executables

$(context).runs        := $(patsubst %.cpp,%.run,$($(context).sources))

run: $(context).run
$(context).run: $($(context).runs)
$($(context).runs): %.run : %.$(ending)
	./$< 

$(context).silent_runs := $(patsubst %.cpp,%.silent_run,$($(context).sources))

silent_run: $(context).silent_run
$(context).silent_run: $($(context).silent_runs)
$($(context).silent_runs): %.silent_run : %.$(ending)
	./$< > /dev/null 

PHONY: $(context).run $(context).silent_run $($(context).runs) $($(context).silent_runs)

# 2> /dev/null