

default:
	cd ./combinatorics && make && cd ..
	cd ./operators && make && cd ..
	cd ./dense && make && cd ..
	cd ./sparse && make && cd ..
	cd ./solver && make && cd ..
	cd ./mesh && make && cd ..
	cd ./vtk && make && cd ..
	cd ./tests && make && cd ..
	echo "finished building." 
    
clean:
	cd ./combinatorics && make clean && cd ..
	cd ./operators && make clean && cd ..
	cd ./dense && make clean && cd ..
	cd ./sparse && make clean && cd ..
	cd ./solver && make clean && cd ..
	cd ./mesh && make clean && cd ..
	cd ./vtk && make clean && cd ..
	cd ./tests && make clean && cd ..
	make -f makefile.clean clean
	echo "finished cleaning." 
    
