#!/usr/bin/env python

import math, urllib, sys, os, traceback
import rpy2.robjects as ro
from metagenome import Metagenome
from ipyTools import *
from collections import defaultdict

class AnalysisSet(object):
    def __init__(self, ids=[], auth=None, adir=None, def_name=None):
        if adir is None:
            adir = random_str()
        self._dir  = adir
        self._path = Ipy.NB_DIR+'/'+adir
        self._auth = auth
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name
        for tax in Ipy.TAX_SET:
            values = {}
            for val in Ipy.VALUES:
                values[val] = self._get_analysis(ids, 'organism', tax, val, 'M5NR')
            setattr(self, tax, values)
        for ont in Ipy.ONT_SET:
            values = {}
            for val in Ipy.VALUES:
                values[val] = self._get_analysis(ids, 'function', ont, val, 'Subsystems')
            setattr(self, ont, values)

    def dump(self):
        for tax in Ipy.TAX_SET:
            tax_set = self.__getitem__(tax)
            for analysis in tax_set.itervalues():
                fname = self._path+'/'+analysis.id+'.biom'
                analysis.dump(fname, fformat='biom')
        for ont in Ipy.ONT_SET:
            ont_set = self.__getitem__(ont)
            for analysis in ont_set.itervalues():
                fname = self._path+'/'+analysis.id+'.biom'
                analysis.dump(fname, fformat='biom')

    def _get_analysis(self, ids, annotation, level, result_type, source):
        # this needs to be created same way as matrix api builds it
        matrix_id = "_".join(sorted(ids))+"_"+"_".join([annotation, level, source, result_type])
        matrix_id += "_%d_%d_%d"%(Ipy.MATRIX['e_val'], Ipy.MATRIX['ident'], Ipy.MATRIX['alen'])
        # load from client cache if exists
        biom_file = self._path+'/'+matrix_id+'.biom'
        if os.path.isfile(biom_file):
            return Analysis(bfile=biom_file, auth=self._auth)
        else:
            keyArgs = Ipy.MATRIX
            keyArgs['ids'] = ids
            keyArgs['annotation'] = annotation
            keyArgs['level'] = level
            keyArgs['result_type'] = result_type
            keyArgs['source'] = source
            if self._auth:
                keyArgs['auth'] = self._auth
            return Analysis(**keyArgs)
    
    def plot_taxon(self, normalize=1, level='domain', parent=None, width=800, height=800, title="", legend=True):
        children = get_taxonomy(level, parent) if parent is not None else None
        keyArgs = { 'normalize': normalize,
                    'ptype': 'row',
                    'width': width,
                    'height': height,
                    'x_rotate': '0',
                    'title': title,
                    'legend': legend,
                    'subset': children }
        if child_level(level, htype='taxonomy'):
            click_opts = (self.defined_name, child_level(level, htype='taxonomy'), normalize, width, height, title, 'True' if legend else 'False')
            keyArgs['onclick'] = "'%s.plot_taxon(level=\"%s\", parent=\"'+params['series']+'\", normalize=%d, width=%d, height=%d, title=\"%s\", legend=%s)'"%click_opts
        self[level]['abundance'].plot_annotation(**keyArgs)
        
    def plot_function(self, normalize=1, level='level1', parent=None, width=800, height=800, title="", legend=True):
        children = get_ontology(level, parent) if parent is not None else None
        keyArgs = { 'normalize': normalize,
                    'ptype': 'row',
                    'width': width,
                    'height': height,
                    'x_rotate': '0',
                    'title': title,
                    'legend': legend,
                    'subset': children }
        if child_level(level, htype='ontology'):
            click_opts = (self.defined_name, child_level(level, htype='ontology'), normalize, width, height, title, 'True' if legend else 'False')
            keyArgs['onclick'] = "'%s.plot_function(level=\"%s\", parent=\"'+params['series']+'\", normalize=%d, width=%d, height=%d, title=\"%s\", legend=%s)'"%click_opts
        self[level]['abundance'].plot_annotation(**keyArgs)

class Analysis(object):
    def __init__(self, ids=[], annotation=None, level=None, result_type=None, source=None, e_val=None, ident=None, alen=None, filters=[], filter_source=None, biom=None, bfile=None, auth=None):
        self._auth = auth
        if (biom is None) and (bfile is None):
            self.biom = self._get_matrix(ids, annotation, level, result_type, source, e_val, ident, alen, filters, filter_source)
        elif biom and isinstance(biom, dict):
            self.biom = biom
        elif bfile and os.path.isfile(bfile):
            try:
                bhdl = open(bfile, 'rU')
                self.biom = json.load(bhdl)
                bhdl.close()
            except:
                self.biom = None
        else:
            self.biom = None
        self._init_matrix()
    
    def _init_matrix(self):
        self.id = self.biom['id'] if self.biom else ""
        self.matrix = self.biom['data'] if self.biom else []
        self.numIDs = self.biom['shape'][1] if self.biom else 0
        self.numAnnotations = self.biom['shape'][0] if self.biom else 0
        self.Dmatrix = self._dense_matrix()
        self.Rmatrix = pyMatrix_to_rMatrix(self.Dmatrix, self.numAnnotations, self.numIDs)
        #self.Nmatrix = self._normalize_matrix()
        self.Nmatrix = self.Rmatrix
        self.alpha_diversity = None
        self.rarefaction     = None
    
    def _get_matrix(self, ids, annotation, level, result_type, source, e_val, ident, alen, filters, filter_source):
        params = map(lambda x: ('id', x), ids)
        if not annotation:
            annotation = Ipy.MATRIX['annotation']
        if level:
            params.append(('group_level', level))
        if result_type:
            params.append(('result_type', result_type))
        if source:
            params.append(('source', source))
        if e_val:
            params.append(('evalue', str(e_val)))
        if ident:
            params.append(('identity', str(ident)))
        if alen:
            params.append(('length', str(alen)))
        if len(filters) > 0:
            params.extend( map(lambda x: ('filter', x), filters) )
            if filter_source:
                params.append(('filter_source', filter_source))
        if self._auth:
            params.append(('auth', self._auth))
        return obj_from_url( Ipy.API_URL+'matrix/'+annotation+'?'+urllib.urlencode(params, True) )
    
    def dump(self, fname, fformat='biom'):
        if self.biom:
            fhdl = open(fname, 'w')
            if fformat == 'biom':
                json.dump(self.biom, fhdl)
            else:
                ann = self.annotations()
                fhdl.write("\t%s\n"%"\t".join(self.ids()))
                for i, row in enumerate(self.Dmatrix):
                    fhdl.write("%s\t%s\n"%(ann[i], "\t".join(row)))
            fhdl.close()
    
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
        mg = Metagenome(aid, auth=self._auth)
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
    
    def rarefaction(self):
        if self.table_type != 'taxon':
            return None
        if self.rarefaction is None:
            rareFact = defaultdict(list)
            for i, aID in enumerate(self.ids()):
                mg = self.get_id_object(aID)
                if ('rarefaction' in mg.stats) and (len(mg.stats['rarefaction']) > 0):
                    rareFact[aID] = mg.stats['rarefaction']
                    continue
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
        fname  = 'images/boxplot_'+random_str()+'.svg'
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
        fname  = 'images/pco_'+random_str()+'.svg'
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
        fname  = 'images/heatmap_'+random_str()+'.svg'
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
    
    def plot_annotation(self, normalize=1, ptype='row', width=800, height=800, x_rotate='0', title="", legend=True, subset=None, onclick=None):
        names = self.names()
        if not (names and self.Dmatrix):
            return None
        labels = []
        matrix = rMatrix_to_pyMatrix(self.Nmatrix, self.numAnnotations, self.numIDs) if normalize else self.Dmatrix
        colors = google_palette(len(names))
        data   = []
        for i, n in enumerate(names):
            data.append({'name': n, 'data': [], 'fill': colors[i]})
        for r, row in enumerate(matrix):
            if (subset is not None) and (self.biom['rows'][r]['id'] not in subset):
                continue
            labels.append(self.biom['rows'][r]['id'])
            for c, val in enumerate(row):
                data[c]['data'].append(val)
        
        keyArgs = { 'btype': ptype,
                    'width': width,
                    'height': height,
                    'x_labels': json.dumps(labels),
                    'x_labels_rotation': x_rotate,
                    'title': title,
                    'target': random_str(),
                    'show_legend': legend,
                    'legend_position': 'right',
                    'data': data,
                    'onclick': onclick }
        try:
            Ipy.RETINA.graph(**keyArgs)
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
