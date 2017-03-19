

default:
	cd ./combinatorics; make; cd ..
	cd ./operators; make; cd ..
	cd ./dense; make; cd ..
	cd ./sparse; make; cd ..
	cd ./solver; make; cd ..
	cd ./mesh; make; cd ..
	echo "finished building." 
    
clean:
	cd ./combinatorics; make clean; cd ..
	cd ./operators; make clean; cd ..
	cd ./dense; make clean; cd ..
	cd ./sparse; make clean; cd ..
	cd ./solver; make clean; cd ..
	cd ./mesh; make clean; cd ..
	echo "finished cleaning." 
    