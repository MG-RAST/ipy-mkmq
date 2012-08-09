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

    def _build_url(self):
        module = 'matrix/'+self.annType if self.idType == 'metagenome' else 'genome_matrix'
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

    def get_pco(self, method='bray-curtis'):
        if not self.biom:
            return None
        dMatrix = self.dense_matrix()
        rMatrix = pyMatrix_to_rMatrix(dMatrix, self.numAnnotations, self.numIDs)
        keyArgs = {'method': method}
        return ro.r.mpco(rMatrix, **keyArgs)

    def plot_pco(self, filename='tmpPCO.png', pco=None, labels=None):
        if pco is None:
            pco = self.get_pco()
        if labels is None:
            labels = self.ids()
        keyArgs = {'fname': filename, 'labels': ro.r.c(labels)}
        ro.r.render(pco, **keyArgs)
        return filename

    def dense_matrix(self):
        if not self.biom:
            return []
        if self.biom['matrix_type'] == 'dense':
            return self.matrix
        else:
            return sparse_to_dense(self.matrix, self.numAnnotations, self.numIDs)
