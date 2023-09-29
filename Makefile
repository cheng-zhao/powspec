all:
	python setup.py build_ext --inplace
clean:
	rm *.so src/pypowspec.c