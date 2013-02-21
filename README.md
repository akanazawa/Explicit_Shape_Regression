Explicit Shape Regression
====================

Implements the cascaded regressor for face alignment described in the
paper: Cao, Xudong, et al. "Face Alignment by Explicit Shape
Regression." Computer Vision and Pattern Recognition (CVPR), 2012 IEEE
Conference on. IEEE, 2012.

How to run
--------------------

The training images are assumed to all come in same size. These images
should be the output of a face detector. Specify parameters and dataset information
in the config.m file such as number of fiducial points, index of left
eye, right eye etc.

Store the training and test images along with their ground truth
fiducial location in a mat file in the format:
* Is: nRow x nCol x nImg matrix 
* Sgts: 2 x nPoints x nImg matrix

Then run

>do_esr('config')


