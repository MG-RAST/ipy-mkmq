#!/usr/bin/env python

from time import localtime, strftime
from collections import defaultdict
import sys, urllib, urllib2, json
import string, random
import rpy2.robjects as ro

API_URL = 'http://api.metagenomics.anl.gov/api2.cgi/'
COLORS  = [ "#3366cc",
            "#dc3912",
            "#ff9900",
            "#109618",
            "#990099",
            "#0099c6",
            "#dd4477",
            "#66aa00",
            "#b82e2e",
            "#316395",
            "#994499",
            "#22aa99",
            "#aaaa11",
            "#6633cc",
            "#e67300",
            "#8b0707",
            "#651067",
            "#329262",
            "#5574a6",
            "#3b3eac",
            "#b77322",
            "#16d620",
            "#b91383",
            "#f4359e",
            "#9c5935",
            "#a9c413",
            "#2a778d",
            "#668d1c",
            "#bea413",
            "#0c5922",
            "#743411" ]

def google_palette(num):
    if not num:
        return COLORS
    num_colors = []
    for i in range(num):
        c_index = i % len(COLORS);
        num_colors.append( COLORS[c_index] )
    return num_colors

def obj_from_url(url):
    try:
        req = urllib2.Request(url, headers={'Accept': 'application/json'})
        res = urllib2.urlopen(req)
    except urllib2.HTTPError, error:
        print "ERROR (%s): %s"%(url, error.read())
        return None
    if not res:
        print "ERROR (%s): no results returned"%url
        return None
    obj = json.loads(res.read())
    if not obj:
        print "ERROR (%s): return structure not valid json format"%url
        return None
    return obj

def slice_column(matrix, index):
    data = []
    for row in matrix:
        data.append(row[index])
    return data

def sparse_to_dense(sMatrix, rmax, cmax):
    dMatrix = [[0 for i in range(cmax)] for j in range(rmax)]
    for sd in sMatrix:
        r, c, v = sd
        dMatrix[r][c] = v
    return dMatrix

def pyMatrix_to_rMatrix(matrix, rmax, cmax):
    if len(matrix) == 0:
        return None
    mList = []
    for i in range(cmax):
        cList = map(lambda x: x[i], matrix)
        mList.extend(cList)
    return ro.r.matrix(ro.IntVector(mList), nrow=rmax)

def random_str(size=8):
    chars = string.ascii_letters + string.digits
    return ''.join(random.choice(chars) for x in range(size))

def merge_biom(b1, b2):
    if b1 and b2 and (b1['type'] == b2['type']) and (b1['matrix_type'] == b2['matrix_type']) and (b1['matrix_element_type'] == b2['matrix_element_type']) and (b1['matrix_element_value'] == b2['matrix_element_value']):
        mBiom = { "generated_by": b1['generated_by'],
                   "matrix_type": 'dense',
                   "date": strftime("%Y-%m-%dT%H:%M:%S", localtime()),
                   "data": [],
                   "rows": [],
                   "matrix_element_value": b1['matrix_element_value'],
                   "matrix_element_type": b1['matrix_element_type'],
                   "format_url": "http://biom-format.org",
                   "format": "Biological Observation Matrix 1.0",
                   "columns": [],
                   "id": b1['id']+'_'+b2['id'],
                   "type": b1['type'],
                   "shape": [] }
        cols, rows = merge_matrix_info(b1['columns'], b2['columns'], b1['rows'], b2['rows'])
        merge_func = merge_sparse if b1['type'] == 'sparse' else merge_dense
        mCol, mRow, mData = merge_func([b1['data'], b2['data']], cols, rows)
        mBiom['columns'] = mCol
        mBiom['rows'] = mRow
        mBiom['data'] = mData
        mBiom['shape'] = [ len(mRow), len(mCol) ]
        return mBiom
    else:
        sys.stderr.write("The inputed biom objects are not compatable for merging\n")
        return None

def merge_matrix_info(c1, c2, r1, r2):
    ## merge columns, skip duplicate
    cm = {}
    for i, c in enumerate(c1):
        cm[ c[id] ] = [0, i, c]
    for i, c in enumerate(c2):
        if c[id] in cm:
            continue
        cm[ c[id] ] = [1, i, c]
    ## merge rows
    rm = defaultdict(list)
    for i, r in enumerate(r1):
        rm[ r[id] ].append( [0, i, r] )
    for i, r in enumerate(r2):
        rm[ r[id] ].append( [1, i, r] )
    return cm.values(), rm.values()

def merge_sparse(data, cols, rows):
    for i in range(len(data)):
        data[i] = sparse_to_dense(data[i], len(cols), len(rows))
    return merge_dense(data, cols, rows)
    
def merge_dense(data, cols, rows):
    cm = map(lambda x: x[2], cols)
    rm = map(lambda x: x[0][2], rows)
    mm = [[0 for i in range(len(cols))] for j in range(len(rows))]
    for i, rset in enumerate(rows):
        for r in rset:
            for j, c in enumerate(cols):
                if r[0] == c[0]:
                    mm[i][j] += data[ r[0] ][ r[1] ][ c[1] ]
    return cm, rm, mm

def get_taxonomy(level='species', parent=None):
    params = []
    if parent is not None:
        params.append(('parent_name', parent))
        params.append(('parent_level', level))
    else:
        params.append(('min_level', level))
    return obj_from_url(API_URL+'m5nr/taxonomy?'+urllib.urlencode(params, True))

