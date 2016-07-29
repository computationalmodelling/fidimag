for ((i = 1; i <= 16; i++)); do
	export OMP_NUM_THREADS=$i
        echo $OMP_NUM_THREADS
	python test_cvode_serial.py
	python test_cvode_parallel.py
done
