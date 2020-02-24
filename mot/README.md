# Multi-Object Tracking

A nice intro to the subject can be found here:
https://www.learnopencv.com/object-tracking-using-opencv-cpp-python/

The two key steps are *Detection* and *Tracking*.  This is a very active
field of research and there are many implementations of each.  We chose
two software projects to represent aspects of the problem:

## GoogLeNet-Inception

https://research.google/pubs/pub43022/
https://github.com/conan7882/GoogLeNet-Inception

This is an implementation of an image classification algorithm called
Inception, described in the above paper.  It applies a pre-made model
to photos and outputs scores for how likely the photo is to be one of
various objects (cat, panda, vehicle, etc).

# Deep SORT

https://arxiv.org/abs/1703.07402
https://github.com/nwojke/deep_sort

This is an implementation of the tracking algorithm described in the
above paper.  It takes detection data and image data as input, and
attempts to track the detections over time and space.

NOTE: If you edit the Makefile and change --display=False to True,
the app will pop up a movie window and draw boxes around pedestrians
as it runs.
