clean:
	-rm -f *.o
	-rm -f $(TARGET)

clean_all:
	make clean
	make pyclean
	make matclean

pyclean:
	-rm -f *.so
	-rm -rf *.egg-info*
	-rm -rf ./tmp/
	-rm -rf ./build/

matclean:
	-rm -rf ./matlabbuild/

python:
	sudo make pyclean
	sudo pip install -e ../cNetworkDiff --no-binary :all:
