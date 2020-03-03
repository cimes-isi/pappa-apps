run: run-detect run-track

# run things
run-detect: venv/stamp GoogLeNet-Inception/data/dropbox/googlenet.npy
	. venv/bin/activate && cd GoogLeNet-Inception/examples && time python3 inception_pretrained.py

run-track: venv/stamp mot16/train/MOT16-04/seqinfo.ini deep_sort/resources/detections/MOT16_train/MOT16-04.npy
	# change --display to True to see the tracking happen visually.
	. venv/bin/activate && cd deep_sort && time python3 deep_sort_app.py \
	                                       --sequence_dir=../mot16/train/MOT16-04 \
	                                       --detection_file=./resources/detections/MOT16_train/MOT16-04.npy \
	                                       --min_confidence=0.3 \
	                                       --nn_budget=100 \
	                                       --display=False


# get googlenet data
GoogLeNet-Inception/data/dropbox/googlenet.zip:
	mkdir -p GoogLeNet-Inception/data/dropbox
	wget -c --max-redirect=20 "https://www.dropbox.com/sh/axnbpd1oe92aoyd/AADpmuFIJTtxS7zkL_LZrROLa?dl=1" -O GoogLeNet-Inception/data/dropbox/googlenet.zip.tmp
	mv GoogLeNet-Inception/data/dropbox/googlenet.zip.tmp GoogLeNet-Inception/data/dropbox/googlenet.zip

GoogLeNet-Inception/data/dropbox/googlenet.npy: GoogLeNet-Inception/data/dropbox/googlenet.zip
	cd GoogLeNet-Inception/data/dropbox && unzip googlenet.zip || true
	touch GoogLeNet-Inception/data/dropbox/googlenet.npy


# get deep_sort data
mot16/MOT16.zip:
	mkdir -p mot16
	wget -c https://motchallenge.net/data/MOT16.zip -O mot16/MOT16.zip.tmp
	mv mot16/MOT16.zip.tmp mot16/MOT16.zip

deep_sort/resources/networks/mars-small128.pb:
	mkdir -p deep_sort/resources/networks
	wget 'https://docs.google.com/uc?export=download&id=1bB66hP9voDXuoBoaCcKYY7a8IYzMMs4P' -O deep_sort/resources/networks/mars-small128.pb.tmp
	mv deep_sort/resources/networks/mars-small128.pb.tmp deep_sort/resources/networks/mars-small128.pb

mot16/train/MOT16-04/seqinfo.ini: mot16/MOT16.zip
	cd mot16 && unzip MOT16.zip
	touch mot16/train/MOT16-04/seqinfo.ini

deep_sort/resources/detections/MOT16_train/MOT16-04.npy: venv/stamp mot16/train/MOT16-04/seqinfo.ini deep_sort/resources/networks/mars-small128.pb
	mkdir -p deep_sort/resources/detections
	. venv/bin/activate && cd deep_sort && python3 tools/generate_detections.py \
	    --model=resources/networks/mars-small128.pb \
	    --mot_dir=../mot16/train \
	    --output_dir=./resources/detections/MOT16_train


# set up python virtual environment
venv/stamp:
	python3 -m virtualenv -p python3 venv
	. venv/bin/activate && pip3 install -r requirements.txt
	touch venv/stamp
