#!/usr/bin/env python

import math, urllib
import rpy2.robjects as ro
import metagenome
from ipyTools import *
from collections import defaultdict

class Analysis:
    def __init__(self, ids=[], annotation='organism', level=None, resultType=None, source=None, biom=None):
        ## load matR
        ro.r('library(matR)')
        if biom is None:
            self.biom = self._get_matrix(ids, annotation, level, resultType, source)
        else:
            self.biom = biom
        self._init_matrix()
    
    def _init_matrix(self):
        self.matrix = self.biom['data'] if self.biom else []
        self.numIDs = self.biom['shape'][1] if self.biom else 0
        self.numAnnotations = self.biom['shape'][0] if self.biom else 0
        self.Dmatrix = self._dense_matrix()
        self.Rmatrix = pyMatrix_to_rMatrix(self.Dmatrix, self.numAnnotations, self.numIDs)
        self.Nmatrix = self._normalize_matrix()
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
        return obj_from_url( API_URL+'matrix/'+annotation+'?'+urllib.urlencode(params, True) )
    
    def ids(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['columns'])
    
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
        if not self.biom:
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
        mg = metagenome.Metagenome(aid)
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
    
    def get_pco(self, method='bray-curtis', normalize=1):
        keyArgs = {'method': method}
        matrix  = self.Nmatrix if normalize else self.Rmatrix
        if not matrix:
            return None
        return ro.r.mpco(matrix, **keyArgs)
    
    def plot_pco(self, filename=None, pco=None, labels=None, title=None):
        if pco is None:
            pco = self.get_pco()
        if labels is None:
            labels = self.ids()
        if filename is None:
            filename = 'pco_'+random_str()+'.png'
        if title is None:
            title = 'PCoA'
        keyArgs = {'fname': filename, 'labels':ro.StrVector(labels), 'main': title}
        ro.r.render(pco, **keyArgs)
        return filename
    
    def plot_heatmap(self, filename=None, normalize=1, labels=None, title=None):
        matrix = self.Nmatrix if normalize else self.Rmatrix
        if not matrix:
            return None
        if labels is None:
            labels = self.ids()
        if filename is None:
            filename = 'heatmap_'+random_str()+'.jpeg'
        if title is None:
            title = 'HeatMap'
        keyArgs = {'image_out': filename, 'labCol': ro.StrVector(labels), 'image_title': title, 'col_lab_mult': 1.2, 'margins': ro.r.c(9,1)}
        ro.r.mheatmap(matrix, **keyArgs)
        return filename
    
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
