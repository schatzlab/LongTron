#!/bin/env python
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import svm
from sklearn.cross_validation import cross_val_score
from numpy import genfromtxt

#cut -f 1,11,12,14,18,19,21- all.tfeatures | tail -n+2 | perl -ne 'BEGIN { %h=("non-recurrent"=>0,"novel"=>1,"problem-free"=>2,"recurrent"=>3); } chomp; @f=split(/\t/,$_); $n=shift(@f); @a=(0,0,0,0); $a[$h{$n}]=1; print "$n\t".join("\t",@f)."\n"; print STDERR "$n\t".join("\t",@a)."\n";' > all.tfeatures.all_x_vectors 2> all.tfeatures.all_y_vectors
if __name__ == '__main__':
    x = genfromtxt("training.x",delimiter='\t')
    y = genfromtxt("training.y",delimiter='\t')

    #only recurrent vs. problem-free
    x2 = genfromtxt("training.2.x",delimiter='\t')
    y2 = genfromtxt("training.2.y",delimiter='\t')
    
    #no novels
    x3 = genfromtxt("training.3.x",delimiter='\t')
    y3 = genfromtxt("training.3.y",delimiter='\t')

    #binary break down
    yb = genfromtxt("binary.y.training",delimiter='\t')
    #permuted labels
    yr = genfromtxt("training.y.shuf",delimiter='\t')
    
    xv = genfromtxt("validation.x",delimiter='\t')
    yv = genfromtxt("validation.y",delimiter='\t')
    
    xv2 = genfromtxt("validation.2.x",delimiter='\t')
    yv2 = genfromtxt("validation.2.y",delimiter='\t')
    
    xv3 = genfromtxt("validation.3.x",delimiter='\t')
    yv3 = genfromtxt("validation.3.y",delimiter='\t')
    
    ybv = genfromtxt("binary.y.validation",delimiter='\t')
    
    xt = genfromtxt("testing.x",delimiter='\t')
    yt= genfromtxt("testing.y",delimiter='\t')

    crf = RandomForestClassifier(n_estimators=100, n_jobs=4)
    c = crf.fit(x,y)
    c.score(xv,yv)
    
    crf = ExtraTreesClassifier(n_estimators=100, n_jobs=4)
    c = crf.fit(x,y)
    c.score(xv,yv)

    cvscores = cross_val_score(crf, x, y, cv=10)

    srf = svm.SVC(verbose=True)
    s = srf.fit(x,y)
    s.predict(xv,yv)
