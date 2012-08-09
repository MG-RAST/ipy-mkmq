#!/usr/bin/env python

import urllib, urllib2, json
import rpy2.robjects as ro

API_URL = 'http://api.metagenomics.anl.gov/'

def obj_from_url(url):
    try:
        req = urllib2.Request(url, headers={'Accept': 'application/json'})
        res = urllib2.urlopen(req)
    except urllib2.HTTPError, error:
        print "ERROR: ", error.read()
        return None
    if not res:
        print "ERROR: no results returned"
        return None
    obj = json.loads(res.read())
    if not obj:
        print "ERROR: return structure not valid json format"
        return None
    return obj

def slice_column(matrix, index):
    data = []
    for row in matrix:
        data.append(row[index])
    return data

def sparse_to_dense(sMatrix, rmax, cmax):
    dMatrix = []
    for r in range(rmax):
        cols = list([0 for x in range(cmax)])
        dMatrix.append(cols)
    for sd in sMatrix:
        r, c, v = sd
        dMatrix[r][c] = v
    return dMatrix

def pyMatrix_to_rMatrix(matrix, rmax, cmax):
    mList = []
    for i in range(cmax):
        cList = map(lambda x: x[i], matrix)
        mList.extend(cList)
    return ro.r.matrix(ro.IntVector(mList), nrow=rmax)
