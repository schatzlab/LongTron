#!/bin/env python
import sys
import os
import re

import cPickle as pickle

from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt
from sklearn.externals import joblib

import matplotlib.pyplot as mp
from sklearn.metrics import roc_curve, auc

fname = sys.argv[1]
#training.bam.bed.rl.nX3.minX2.mq.rm.ot.ed.td.logsX3.sm.gc.lsX10
pos = fname.find('rl')
fields_ = fname[pos:].split('.')
fields = []
p = re.compile(r'([^X]+)X(\d+)')
for f in fields_:
    m = p.search(f)
    if m:
        label = m.group(1)
        multiplier = m.group(2)
        fields.extend(int(multiplier)*[label])
    else:
        fields.append(f)

print fields
print len(fields)

x = genfromtxt("training.x",delimiter='\t')
#x = genfromtxt("training.x.1k",delimiter='\t')
#y = genfromtxt("training.y",delimiter='\t')

#binary break down
yb = genfromtxt("binary.y.training",delimiter='\t')
#yb = genfromtxt("binary.y.training.1k",delimiter='\t')
#permuted labels
#yr = genfromtxt("training.y.shuf",delimiter='\t')

xv = genfromtxt("validation.x",delimiter='\t')
#xv = genfromtxt("validation.x.1k",delimiter='\t')
#yv = genfromtxt("validation.y",delimiter='\t')

ybv = genfromtxt("binary.y.validation",delimiter='\t')
#ybv = genfromtxt("binary.y.validation.1k",delimiter='\t')

xt = genfromtxt("testing.x",delimiter='\t')
ybt = genfromtxt("binary.y.testing",delimiter='\t')
#xt = genfromtxt("testing.x.1k",delimiter='\t')
#ybt = genfromtxt("binary.y.testing.1k",delimiter='\t')

if not os.path.exists('./trained_models'):
    os.mkdir('./trained_models')

#crf = RandomForestClassifier(n_estimators=100, n_jobs=8)
#c = crf.fit(x,y)
#joblib.dump(c,"./trained_models/%s.rfc_4_trained" % (fname))
#sys.stdout.write("validation: %s\n" % str(c.score(xv,yv)))
#sys.stdout.write("importance: ")
#[sys.stdout.write(" %s:%s " % (fields[i],str(x_))) for (i,x_) in enumerate(c.feature_importances_)]
#sys.stdout.write("\n")

bcrf = RandomForestClassifier(n_estimators=100, n_jobs=8)
bc = bcrf.fit(x,yb)
joblib.dump(bc,"./trained_models/%s.rfc_2_trained" % (fname))

vscore = bc.score(xv,ybv)
vprobs = bc.predict_proba(xv)
sys.stdout.write("binary validation: %s\n" % str(vscore))

sys.stdout.write("importance: ")
[sys.stdout.write(" %s:%s " % (fields[i],str(x_))) for (i,x_) in enumerate(bc.feature_importances_)]
sys.stdout.write("\n")

fout = open("bvprobs.pkl","wb")
pickle.dump(vprobs, fout)
fout.close()

#os.symlink("./%s.rfc_4_trained" % (fname),'./trained_models/short.rfc_4_trained')
#os.symlink("./%s.rfc_2_trained" % (fname),'./trained_models/short.rfc_2_trained')

tscore = bc.score(xt,ybt)
tprobs = bc.predict_proba(xt)
sys.stdout.write("binary test: %s\n" % str(tscore))
fout = open("btprobs.pkl","wb")
pickle.dump(tprobs, fout)
fout.close()

#ROC plotting
tpr = dict()
fpr = dict()
rauc = dict()
fpr[0], tpr[0], _ = roc_curve(ybt, tprobs[:,1]) 
#fpr[0], tpr[0], _ = roc_curve(ybv, vprobs[:,1])
rauc[0] = auc(fpr[0], tpr[0])
#savefig
mp.figure()
lw = 2
mp.plot(fpr[0], tpr[0], color='darkorange', lw=lw, label='ROC curve (auc = %0.3f)' % rauc[0])
#mp.plot(fpr[0], tpr[0], color='darkorange', lw=lw, label='ROC curve')
mp.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
mp.xlim([0.0, 1.0])
mp.ylim([0.0, 1.05])
mp.xlabel('FPR')
mp.ylabel('TPR')
mp.title('PacBio FL Binary Testing ROC')
mp.legend(loc='lower right')
mp.savefig('pbfl.b.roc.png', bbox_inches='tight')
