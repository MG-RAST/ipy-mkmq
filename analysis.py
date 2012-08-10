#!/usr/bin/env python

import rpy2.robjects as ro
from ipyTools import *

class Analysis:
    def __init__(self, ids, idType='metagenome', annotation='organism', level=None, resultType=None, source=None):
        self.inputIDs   = ids
        self.idType     = idType
        self.annotation = annotation
        self.resultType = resultType
        self.source     = source
        self.level      = level
        self.matrixUrl  = self._build_url()
        self.biom       = obj_from_url(self.matrixUrl)
        self.matrix     = self.biom['data'] if self.biom else []
        self.numIDs     = self.biom['shape'][1] if self.biom else 0
        self.numAnnotations = self.biom['shape'][0] if self.biom else 0
        dMatrix = self.dense_matrix()
        self.Rmatrix  = pyMatrix_to_rMatrix(dMatrix, self.numAnnotations, self.numIDs) if len(dMatrix) > 0 else None
        self.NRmatrix = self.normalize_matrix()

    def _build_url(self):
        module = 'matrix/'+self.annotation if self.idType == 'metagenome' else 'genome_matrix'
        params = [ ('format', 'biom'), ('id', self.inputIDs) ]
        if self.resultType and (self.idType == 'metagenome'):
            params.append(('result_type', self.resultType))
        if self.source and (self.idType == 'metagenome'):
            params.append(('source', self.source))
        if self.level:
            params.append(('group_level', self.level))
        return API_URL+module+'?'+urllib.urlencode(params, True)

    def ids(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['columns'])

    def annotations(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['rows'])

    def data_for_id(self, id):
        if not self.biom:
            return None
        items = self.ids()
        index = items.index(id)
        return slice_column(self.matrix, index)

    def normalize_matrix(self):
        if self.Rmatrix:
            return ro.r.normalize(self.Rmatrix)
        else:
            return None

    def get_pco(self, method='bray-curtis', normalize=1):
        keyArgs = {'method': method}
        matrix  = self.NRmatrix if normalize else self.Rmatrix
        if not matrix:
            return None
        return ro.r.mpco(matrix, **keyArgs)

    def plot_pco(self, filename=None, pco=None, labels=None, title=None):
        if pco is None:
            pco = self.get_pco()
        if labels is None:
            labels = self.ids()
        if filename is None:
            filename = 'images/pco_'+random_str()+'.png'
        if title is None:
            title = 'PCoA'
        keyArgs = {'fname': filename, 'labels':ro.StrVector(labels), 'main': title}
        ro.r.render(pco, **keyArgs)
        return filename

    def plot_heatmap(self, filename=None, normalize=1, labels=None, title=None):
	matrix = self.NRmatrix if normalize else self.Rmatrix
        if not matrix:
            return None
        if labels is None:
            labels = self.ids()
        if filename is None:
            filename = 'images/heatmap_'+random_str()+'.jpeg'
        if title is None:
            title = 'HeatMap'
        keyArgs = {'image_out': filename, 'labCol': ro.StrVector(labels), 'image_title': title, 'col_lab_mult': 1.2, 'margins': ro.r.c(9,1)}
        ro.r.mheatmap(matrix, **keyArgs)
	return filename

    def dense_matrix(self):
        if not self.biom:
            return []
        if self.biom['matrix_type'] == 'dense':
            return self.matrix
        else:
            return sparse_to_dense(self.matrix, self.numAnnotations, self.numIDs)
