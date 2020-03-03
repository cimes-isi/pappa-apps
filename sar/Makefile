run: venv/stamp
	. venv/bin/activate && cd original/examples && time python3 AutoFocus_demo.py --batch

# set up python virtual environment
venv/stamp:
	python3 -m virtualenv -p python3 venv
	. venv/bin/activate && pip3 install -r requirements.txt
	touch venv/stamp
