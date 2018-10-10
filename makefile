

default:
	cd ./combinatorics && $(MAKE) 
	cd ./operators && $(MAKE) 
	cd ./dense && $(MAKE) 
	cd ./sparse && $(MAKE) 
	cd ./solver && $(MAKE) 
	cd ./mesh && $(MAKE) 
	cd ./vtk && $(MAKE)
	cd ./tests && $(MAKE)
	echo "finished building." 
    
clean:
	cd ./combinatorics && $(MAKE) clean
	cd ./operators && $(MAKE) clean
	cd ./dense && $(MAKE) clean
	cd ./sparse && $(MAKE) clean
	cd ./solver && $(MAKE) clean
	cd ./mesh && $(MAKE) clean
	cd ./vtk && $(MAKE) clean
	cd ./tests && $(MAKE) clean
	$(MAKE) -f $(MAKE)file.clean clean
	@echo "finished cleaning." 
    
