

$(context).sources     := $(wildcard $(contextdir)/*cpp)
$(context).executables := $(patsubst %.cpp,%.out,$($(context).sources))
$(context).runs        := $(patsubst %.cpp,%.run,$($(context).sources))

# $(context).runs        := $(patsubst %.out,%.run,$($(context).executables))


# .PHONY: %.run

$($(context).runs): %.run : %.out
#	./$< > /dev/null
	./$< 

$(context).run: $($(context).runs)

run: $(context).run
