# # Example usage:
# # 
# # 
# # include ../../common.recipe.mk 
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



$(context).sources := $(sort $(wildcard $(contextdir)/*.cpp))

$(context).outs    := $(patsubst %.cpp,%.out,$($(context).sources))

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

$($(context).depdir):
	@mkdir -p $@






$($(context).outs): mycontext    := $(context)
$($(context).outs): mycontextdir := $(contextdir)
$($(context).outs): $(contextdir)/%.out: $(contextdir)/%.cpp | $($(context).depdir)
# 	@ echo $(mycontext)
# 	@ echo $(mycontextdir)
# 	@ echo target $@
# 	@ echo target $<
# 	@ echo target $^
# 	@ echo contextdir $(mycontextdir)
# 	@ echo depdir $($(mycontext).depdir)
# 	@ echo $($(mycontext).include)
# 	@ echo $($(mycontext).rpath)
# 	@ echo $($(mycontext).lib)
	@g++ -MM $(mycontextdir)/$*.cpp -MT $@ -MF $($(mycontext).depdir)/$*.d
	$(CXX) $(CXXFLAGS_EXECUTABLE) $(CPPFLAGS) $< $($(mycontext).include) $($(mycontext).rpath) $($(mycontext).mylib) -o $@ $(LDLIBS)


-include $($(context).dependencies)

$($(context).outs): $(contextdir)/%.out: $(contextdir)/makefile 
$($(context).outs): $(contextdir)/%.out: $(projectdir)/makefile 
$($(context).outs): $(contextdir)/%.out: $(projectdir)/*.mk $(projectdir)/tests/*.mk

$(context).tests: $($(context).outs)
