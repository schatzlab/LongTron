#!/bin/env python
import sys
import re
from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt
from sklearn.externals import joblib

#model filename
fname = sys.argv[1]
c = joblib.load('trained_models/%s.rfc_trained4' % fname)
xv = genfromtxt("validation.x",delimiter='\t')
yv = genfromtxt("validation.y",delimiter='\t')
sys.stdout.write("validation: %s\n" % str(c.score(xv,yv)))

bc = joblib.load('trained_models/%s.rfc_trained2' % fname)
ybv = genfromtxt("binary.y.validation",delimiter='\t')
sys.stdout.write("binary validation: %s\n" % str(bc.score(xv,ybv)))
