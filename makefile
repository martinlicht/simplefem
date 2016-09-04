

default:
	cd ./combinatorics; make; cd ..
	cd ./operators; make; cd ..
	cd ./mesh; make; cd ..
	echo "finished building." 
    
clean:
	cd ./combinatorics; make clean; cd ..
	cd ./operators; make clean; cd ..
	cd ./mesh; make clean; cd ..
	echo "finished cleaning." 
    