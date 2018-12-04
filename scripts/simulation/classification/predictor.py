#!/bin/env python
import sys
import re
from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt
from sklearn.externals import joblib
from io import StringIO

MAX_BUF=1000

#previously trained model filename
model_fname = sys.argv[1]
model_class_num = int(sys.argv[2])
#features of real data to predict classes for
#filename e.g.: SRR1163655.sorted.bam.bed.rl.nX3.minX2.mq.rm.sr.snps.ot.gc.umap.ed.td.logsX3.sm.sdX2.lmX4.lsX10
#real_data_fname = sys.argv[2]
c = joblib.load('./trained_models/%s.rfc_%d_trained' % (model_fname, model_class_num))
#buffer up a bunch of lines, then do the prediction
buf = ""
i = 0
for line in sys.stdin:
    if i == MAX_BUF:
        feature_fields = genfromtxt(StringIO(buf.decode('utf-8')), delimiter='\t')
        xvy = c.predict(feature_fields)
        [sys.stdout.write('\t'.join([str(int(z)) for z in y])+'\n') for y in xvy.tolist()]
        i = 0
        buf = ""
    buf += line
    i += 1
if i > 0:
    feature_fields = genfromtxt(StringIO(buf.decode('utf-8')), delimiter='\t')
    xvy = c.predict(feature_fields)
    [sys.stdout.write('\t'.join([str(int(z)) for z in y])+'\n') for y in xvy.tolist()]
    
sys.stdout.write("%d-class prediction finished\n" % (model_class_num))
