run: venv/stamp
	. venv/bin/activate && cd original/examples && time python3 AutoFocus_demo.py --batch

run-bp: venv/stamp
	. venv/bin/activate && cd original/examples && time python3 BP_demo.py --batch

run-ffbp: venv/stamp
	. venv/bin/activate && cd original/examples && time python3 FFBP_demo.py --batch

run-ffbpmp: venv/stamp
	. venv/bin/activate && cd original/examples && time python3 FFBPmp_demo.py --batch

run-ffbpmpi: venv/stamp
	. venv/bin/activate && cd original/examples && time mpirun -n 4 python3 FFBPmpi_demo.py --batch

run-ffbpmpi-sq: venv/stamp
	. venv/bin/activate && cd original/examples && time mpirun -n 4 python3 FFBPmpi_sq_demo.py --batch

run-ffbpmpi-sq-opt: venv/stamp
	. venv/bin/activate && cd original/examples && time mpirun -n 4 python3 FFBPmpi_sq_opt_demo.py --batch

run-ffbpmpi-row-opt: venv/stamp
	. venv/bin/activate && cd original/examples && time mpirun -n 4 python3 FFBPmpi_row_opt_demo.py --batch

exp-afrl: venv/stamp
	. venv/bin/activate && cd original/examples && python3 AFRL_export.py --batch

exp-dirsig: venv/stamp
	. venv/bin/activate && cd original/examples && python3 DIRSIG_export.py --batch

exp-sandia: venv/stamp
	. venv/bin/activate && cd original/examples && python3 Sandia_export.py --batch

# set up python virtual environment
venv/stamp: requirements.txt
	python3 -m virtualenv -p python3 venv
	. venv/bin/activate && pip3 install -r requirements.txt
	touch venv/stamp
