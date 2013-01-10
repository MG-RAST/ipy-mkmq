#!/usr/bin/env python

import math, urllib, sys, os
import rpy2.robjects as ro
from metagenome import Metagenome
from ipyTools import *
from collections import defaultdict

class Analysis(object):
    def __init__(self, ids=[], annotation='organism', level=None, resultType=None, source=None, biom=None, bfile=None, auth=None):
        self._auth = auth
        if (biom is None) and (bfile is None):
            self.biom = self._get_matrix(ids, annotation, level, resultType, source)
        elif biom and isinstance(biom, dict):
            self.biom = biom
        elif bfile and os.path.isfile(bfile):
            try:
                bhdl = open(bfile, 'rU')
                self.biom = json.load(bhdl)
            except:
                self.biom = None
        else:
            self.biom = None
        self._init_matrix()
    
    def _init_matrix(self):
        self.matrix = self.biom['data'] if self.biom else []
        self.numIDs = self.biom['shape'][1] if self.biom else 0
        self.numAnnotations = self.biom['shape'][0] if self.biom else 0
        self.Dmatrix = self._dense_matrix()
        self.Rmatrix = pyMatrix_to_rMatrix(self.Dmatrix, self.numAnnotations, self.numIDs)
        #self.Nmatrix = self._normalize_matrix()
        self.Nmatrix = self.Rmatrix
        self.alpha_diversity = None
        self.rarefaction     = None
    
    def _get_matrix(self, ids, annotation, level, resultType, source):
        params = map(lambda x: ('id', x), ids)
        if level:
            params.append(('group_level', level))
        if resultType:
            params.append(('result_type', resultType))
        if source:
            params.append(('source', source))
        if self._auth:
            params.append(('auth', self._auth))
        return obj_from_url( API_URL+'matrix/'+annotation+'?'+urllib.urlencode(params, True) )
    
    def ids(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['columns'])
    
    def names(self):
        if not self.biom:
            return None
        return map(lambda x: x['name'], self.biom['columns'])
    
    def annotations(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['rows'])
    
    def table_type(self):
        if self.biom and ('type' in self.biom):
            return self.biom['type'].split(' ')[0].lower()
        else:
            return None
    
    def get_id_data(self, aid):
        if not self.Dmatrix:
            return None
        try:
            items = self.ids()
            index = items.index(aid)
        except (ValueError, AttributeError):
            return None
        return slice_column(self.Dmatrix, index)
    
    def get_id_object(self, aid):
        if not self.biom:
            return None
        try:
            items = self.ids()
            index = items.index(aid)
        except (ValueError, AttributeError):
            return None        
        mg = Metagenome(aid)
        if mg.name is not None:
            return mg
        else:
            return self.biom['columns'][index]
    
    def alpha_diversity(self):
        if self.table_type != 'taxon':
            return None
        if self.alpha_diversity is None:
            alphaDiv = {}
            for i, aID in enumerate(self.ids()):
                col = slice_column(self.Dmatrix, i)
                h1  = 0
                s1  = sum(col)
                if not s1:
                    alphaDiv[aID] = 0
                    continue
                for n in col:
                    p = n/s1
                    if p > 0:
                        h1 += (p * math.log(1/p)) / math.log(2)
                alphaDiv[aID] = 2**h1
            self.alpha_diversity = alphaDiv
        return self.alpha_diversity
    
    def rarefaction(self, id):
        if self.table_type != 'taxon':
            return None
        if self.rarefaction is None:
            rareFact = defaultdict(list)
            for i, aID in enumerate(self.ids()):
                mg = self.get_id_object(id)
                try:
                    nseq = int(mg.stats['sequence_count_raw'])
                    size = int(nseq/1000) if nseq > 1000 else 1
                except (ValueError, KeyError, TypeError, AttributeError):
                    rareFact[aID] = []
                    continue
                nums = slice_column(self.Dmatrix, i)
                lnum = len(nums)
                nums.sort()
                for i in xrange(0, nseq, size):
                    coeff = self._nCr2ln(nseq, i)
                    curr  = 0
                    for n in nums:
                        curr += math.exp(self._nCr2ln(nseq-n, i) - coeff)
                    rareFact[aID].append([i, lnum-curr])
            self.rarefaction = rareFact
        return self.rarefaction
    
    def _nCr2ln(self, n, r):
        c = 1
        if r > n:
            return c
        if (r < 50) and (n < 50):
            for x in xrange(0, r-1):
                c += (c * (n-x)) / (x+1)
            return math.log(c)
        if r <= n:
            c = self._gammaln(n+1) - self._gammaln(r+1) - self._gammaln(n-r)
        else:
            c = -1000
        return c
    
    def _gammaln(self, x):
        if x > 0:
            s = math.log(x)
            return math.log(2 * math.pi) / 2 + x * s + s / 2 - x
        else:
            return 0
    
    def plot_boxplot(self, normalize=1, labels=None, title=''):
        matrix = self.Nmatrix if normalize else self.Rmatrix
        fname  = 'boxplot_'+random_str()+'.svg'
        if not matrix:
            return None
        if labels is None:
            labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        keyArgs = { 'names': ro.StrVector(labels),
                    'main': title,
                    'show.names': True,
                    'las': 2,
                    'outpch': 21,
                    'outcex': 0.5,
                    'cex.lab': 0.8,
                    'boxwex': 0.6,
                    'cex.axis': 0.7 }
        ro.r.svg(fname)
        ro.r.boxplot(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def plot_pco(self, normalize=1, method='bray-curtis', labels=None, title=''):
        matrix = self.Nmatrix if normalize else self.Rmatrix
        fname  = 'pco_'+random_str()+'.svg'
        if not matrix:
            return None
        if labels is None:
            labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        keyArgs = { 'labels': ro.StrVector(labels),
                    'main': title,
                    'method': method,
                    'comp': ro.r.c(1,2,3) }
        ro.r.svg(fname)
        ro.r.pco(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def plot_heatmap(self, normalize=1, labels=None, title=''):
        matrix = self.Nmatrix if normalize else self.Rmatrix
        fname  = 'heatmap_'+random_str()+'.svg'
        if not matrix:
            return None
        if labels is None:
            labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        keyArgs = { 'labCol': ro.StrVector(labels),
                    'labRow': '',
                    'main': title,
                    'cexCol': 0.95,
                    'margins': ro.r.c(8,1) }
        ro.r.svg(fname)
        ro.r.heatmap(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def plot_annotation(self, normalize=1, ptype='column', width=1100, height=400, x_rotate='300', title=None, legend=True):
        labels = self.annotations()
        names  = self.names()
        if not (labels and names and self.Dmatrix):
            return None
        matrix = rMatrix_to_pyMatrix(self.Nmatrix, self.numAnnotations, self.numIDs) if normalize else self.Dmatrix
        colors = google_palette(len(names))
        data   = []
        for i, n in enumerate(names):
            data.append({'name': n, 'data': [], 'fill': colors[i]})
        for row in matrix:
            for i, val in enumerate(row):
                data[i]['data'].append(val)
        keyArgs = { 'btype': ptype,
                    'width': width,
                    'height': height,
                    'x_labels': json.dumps(labels),
                    'x_labels_rotation': x_rotate,
                    'title': title,
                    'target': random_str(),
                    'show_legend': legend,
                    'legend_position': 'right',
                    'data': data }
        try:
            RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            print None
    
    def _normalize_matrix(self):
        if self.Rmatrix:
            return ro.r.normalize(self.Rmatrix)
        else:
            return None
    
    def _dense_matrix(self):
        if not self.biom:
            return []
        if self.biom['matrix_type'] == 'dense':
            return self.matrix
        else:
            return sparse_to_dense(self.matrix, self.numAnnotations, self.numIDs)
