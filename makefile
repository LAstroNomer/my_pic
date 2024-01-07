
res: pic
	./pic

pic: pic.o io initial_conditions first second
	gfortran *.o -o pic

pic.o: pic.f95 io initial_conditions first second
	gfortran -c pic.f95

initial_conditions : initial_conditions.f95
	gfortran -c initial_conditions.f95

io: my_io.f95 initial_conditions
	gfortran -c my_io.f95

first: initial_conditions
	gfortran -c first_step.f95

second: initial_conditions
	gfortran -c second_step.f95

plot:
	python3 plot.py

clean:
	rm -f *.o *.mod *.out *.png
