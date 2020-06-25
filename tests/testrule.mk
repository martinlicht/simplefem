$(context).sources := $(wildcard $(context)/*.cpp)

$(context).outs := $(patsubst %.cpp,%.out,$($(context).sources))

$(context).depdir := $(context)/$(depdir)

$(context).dependencies := $(patsubst $(context)/%.cpp,$(context)/.deps/%.d,$($(context).sources))

$(context).include := $(patsubst %,-L../%,$(affix.$(context)))
linkerprefix=-Wl,
$(context).rpath_t := $(patsubst %,-rpath=$(shell pwd)/../%,$(affix.$(context))) 
$(context).rpath := $(patsubst %,$(linkerprefix)%,$($(context).rpath_t)) 
$(context).link := $(patsubst %,-l%,$(affix.$(context)))

$($(context).depdir): ; @mkdir -p $@

$($(context).outs): $(context)/%.out: $(context)/%.cpp | $($(context).depdir)
	@ #echo $($(context).include)
	@ #echo $($(context).rpath)
	@ #echo $($(context).link)
	@g++ -MM $(context)/$*.cpp -MT $@ -MF $(context)/.deps/$*.d
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $($(context).include) $($(context).rpath) $($(context).link) -o $@ $(LDLIBS)

-include $($(context).dependencies)

$(context).tests: $($(context).outs)
