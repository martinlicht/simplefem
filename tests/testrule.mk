
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



$(context).sources := $(wildcard $(contextdir)/*.cpp)

$(context).outs := $(patsubst %.cpp,%.out,$($(context).sources))

$(context).depdir := $(context)/$(depdir)

$(context).dependencies := $(patsubst $(contextdir)/%.cpp,$(contextdir)/.deps/%.d,$($(context).sources))

$(context).include := $(patsubst %,-L$(projectdir)/%,$(affix.$(context)))
linkerprefix=-Wl,
$(context).rpath_t := $(patsubst %,-rpath=$(pathvar)/%,$(affix.$(context))) 
$(context).rpath := $(patsubst %,$(linkerprefix)%,$($(context).rpath_t)) 
$(context).link := $(patsubst %,-l%,$(affix.$(context)))

$($(context).depdir): ; @mkdir -p $@

$($(context).outs): $(contextdir)/%.out: $(contextdir)/%.cpp | $($(context).depdir)
	@ #echo $($(context).include)
	@ #echo $($(context).rpath)
	@ #echo $($(context).link)
	@g++ -MM $(contextdir)/$*.cpp -MT $@ -MF $(contextdir)/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $($(context).include) $($(context).rpath) $($(context).link) -o $@ $(LDLIBS)

-include $($(context).dependencies)

$(context).tests: $($(context).outs)
