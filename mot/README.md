# Multi-Object Tracking

A nice intro to the subject can be found here:
* https://www.learnopencv.com/object-tracking-using-opencv-cpp-python/

The two key steps are *Detection* and *Tracking*.

Detection is an algorithm which identifies objects within a photo.
One way to do this is by repeatedly running an image classifier on
small parts of the image.  This is rather expensive, and does not
understand time at all.  In some situations, it can be more efficient
to run detection rarely and allow tracking to do most of the work.

Tracking is an algorithm which follows objects over time, across
detections in multiple video frames.  This includes velocity estimation,
and re-identification of objects which were not detected for one or more
frames.  For example, a pedestrian may momentarily pass behind a
telephone pole as they walk along.  Re-identification allows the tracker
to detect whether the two detections are, in fact, the same person.

This is a very active field of research and there are many
implementations of detection and tracking.  We chose two software
projects to represent aspects of the MOT problem:

## GoogLeNet-Inception

* https://research.google/pubs/pub43022/
* https://github.com/conan7882/GoogLeNet-Inception

This is an implementation of an image classification algorithm called
Inception, described in the above paper.  It applies a pre-made model
to photos and outputs scores for how likely the photo is to be one of
various objects (cat, panda, vehicle, etc).

# Deep SORT

* https://arxiv.org/abs/1703.07402
* https://github.com/nwojke/deep_sort

This is an implementation of the tracking algorithm described in the
above paper.  It takes detection data and image data as input, and
attempts to track the detections over time and space.

This project also includes detection code, which is run initially to
provide data input for the tracking step.  The detection step is slower
than the tracking, so it is performed once and reused for many runs of
the tracker.

NOTE: If you edit the Makefile and change --display=False to True,
the app will pop up a movie window and draw boxes around pedestrians
as it runs.
