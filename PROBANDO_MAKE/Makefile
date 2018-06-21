
numbers= ${ARGS} 

all: figura

figura: grafica Factorial.png
	@echo "abriendo figura........"
	open Factorial.png

grafica:  imprime  data.dat 
	@echo "graficando............."
	python grafica.py

imprime: 
	@echo "corriendo programa......"
	for number in $(numbers) ; do \
		python imprime_numeros.py $$number ; \
	done


clean:
	rm *.dat *.png