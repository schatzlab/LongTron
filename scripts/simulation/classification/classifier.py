#!/bin/env python
import sys
import re
from sklearn.ensemble import RandomForestClassifier
from numpy import genfromtxt
from sklearn.externals import joblib

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
y = genfromtxt("training.y",delimiter='\t')

#binary break down
yb = genfromtxt("binary.y.training",delimiter='\t')
#permuted labels
#yr = genfromtxt("training.y.shuf",delimiter='\t')

xv = genfromtxt("validation.x",delimiter='\t')
yv = genfromtxt("validation.y",delimiter='\t')

ybv = genfromtxt("binary.y.validation",delimiter='\t')

crf = RandomForestClassifier(n_estimators=100, n_jobs=8)
c = crf.fit(x,y)
#joblib.dump(c,"%s.rfc_trained4" % (fname))
sys.stdout.write("validation: %s\n" % str(c.score(xv,yv)))
sys.stdout.write("importance: ")
#print c.feature_importances_
[sys.stdout.write(" %s:%s " % (fields[i],str(x_))) for (i,x_) in enumerate(c.feature_importances_)]
sys.stdout.write("\n")

bcrf = RandomForestClassifier(n_estimators=100, n_jobs=8)
bc = bcrf.fit(x,yb)
#joblib.dump(bc,"%s.rfc_trained2" % (fname))
sys.stdout.write("binary validation: %s\n" % str(bc.score(xv,ybv)))
sys.stdout.write("importance: ")
[sys.stdout.write(" %s:%s " % (fields[i],str(x_))) for (i,x_) in enumerate(bc.feature_importances_)]
sys.stdout.write("\n")
