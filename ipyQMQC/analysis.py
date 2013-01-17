#!/usr/bin/env python

import math, urllib, sys, os, traceback
import rpy2.robjects as ro
from metagenome import Metagenome
from ipyTools import *
from collections import defaultdict

class AnalysisSet(object):
    def __init__(self, ids=[], auth=None, cache=None, def_name=None):
        if cache is None:
            cache = random_str()
        self._dir  = cache
        self._path = Ipy.NB_DIR+'/'+cache
        self._auth = auth
        self.all_mgs = ids
        self.display_mgs = self.all_mgs
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

    def set_display_mgs(self, ids=[]):
        if (not ids) or (len(ids) == 0):
            sys.stdout.write("setting %s.display_mgs to all metagenomes in set")%self.defined_name
            self.display_mgs = self.all_mgs
        else:
            self.display_mgs = ids

    def dump(self):
        if not os.path.isdir(self._path):
            os.mkdir(self._path)
        for tax in Ipy.TAX_SET:
            tax_set = getattr(self, tax)
            for analysis in tax_set.itervalues():
                fname = self._path+'/'+analysis.id+'.biom'
                analysis.dump(fname, fformat='biom')
        for ont in Ipy.ONT_SET:
            ont_set = getattr(self, ont)
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
            if Ipy.DEBUG:
                sys.stdout.write("loading %s.biom from cache %s ... \n"%(matrix_id, self._dir))
            return Analysis(bfile=biom_file, auth=self._auth)
        else:
            if Ipy.DEBUG:
                sys.stdout.write("loading %s.biom through api ... \n"%matrix_id)
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
        keyArgs = { 'atype': 'taxonomy',
                    'normalize': normalize,
                    'level': level,
                    'parent': parent,
                    'width': width,
                    'height': height,
                    'title': title,
                    'legend': legend }
        self._plot_annotation(**keyArgs)
        
    def plot_function(self, normalize=1, level='level1', parent=None, width=800, height=800, title="", legend=True):
        keyArgs = { 'atype': 'ontology',
                    'normalize': normalize,
                    'level': level,
                    'parent': parent,
                    'width': width,
                    'height': height,
                    'title': title,
                    'legend': legend }
        self._plot_annotation(**keyArgs)
    
    def _plot_annotation(self, atype='taxonomy', normalize=1, level='domain', parent=None, width=800, height=800, title="", legend=True):
        children = get_hierarchy(htype=atype, level=level, parent=parent)
        keyArgs = { 'normalize': normalize,
                    'ptype': 'row',
                    'width': width,
                    'height': height,
                    'x_rotate': '0',
                    'title': title,
                    'legend': legend,
                    'submg': self.display_mgs,
                    'subset': children }
        if child_level(level, htype=atype):
            click_opts = (self.defined_name, 'taxon' if atype == 'taxonomy' else 'function', child_level(level, htype=atype), normalize, width, height, title, 'True' if legend else 'False')
            keyArgs['onclick'] = "'%s.plot_%s(level=\"%s\", parent=\"'+params['series']+'\", normalize=%d, width=%d, height=%d, title=\"%s\", legend=%s)'"%click_opts
        to_plot = getattr(self, level)
        to_plot['abundance'].plot_annotation(**keyArgs)

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
        self.result_type = self.biom['matrix_element_value'] if self.biom else ""
        self.matrix = self.biom['data'] if self.biom else []
        self.numIDs = self.biom['shape'][1] if self.biom else 0
        self.numAnnotations = self.biom['shape'][0] if self.biom else 0
        self.Dmatrix  = self._dense_matrix()  # count dense matrix
        self.NDmatrix = None  # normalized dense matrix
        self.Rmatrix  = pyMatrix_to_rMatrix(self.Dmatrix, self.numAnnotations, self.numIDs) # R count matrix object
        self.NRmatrix = None  # R normalized matrix object
        if self.result_type == 'abundance':
            self._normalize_matrix() # only normalize abundance counts
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
    
    def dump(self, fname=None, fformat='biom', normalize=0):
        if not fname:
            fname = self.id+'.biom' if fformat == 'biom' else self.id+'.tab'
        if self.biom:
            fhdl = open(fname, 'w')
            if fformat == 'biom':
                json.dump(self.biom, fhdl)
            else:
                matrix = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
                annot  = self.annotations()
                fhdl.write("\t%s\n"%"\t".join(self.ids()))
                for i, row in enumerate(matrix):
                    fhdl.write(annot[i])
                    for r in row:
                        fhdl.write("\t"+str(r))
                    fhdl.write("\n")
            fhdl.close()
        return fname
    
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
    
    def plot_boxplot(self, normalize=1, title=''):
        matrix = self.NRmatrix if normalize and self.NRmatrix else self.Rmatrix
        fname  = Ipy.IMG_DIR+'/boxplot_'+random_str()+'.svg'
        labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        if not matrix:
            return None
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
    
    def plot_pco(self, normalize=1, dist='bray-curtis', title=''):
        matrix = self.NRmatrix if normalize and self.NRmatrix else self.Rmatrix
        fname  = Ipy.IMG_DIR+'/pco_'+random_str()+'.svg'
        labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        if not matrix:
            return None
        keyArgs = { 'labels': ro.StrVector(labels),
                    'main': title,
                    'method': dist,
                    'comp': ro.r.c(1,2,3) }
        ro.r.svg(fname)
        ro.r.pco(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def plot_heatmap(self, normalize=1, title='', source='retina', dist='bray-curtis', clust='ward', width=700, height=600):
        keyArgs = {'normalize': normalize, 'title': title, 'dist': dist, 'clust': clust, 'width': width, 'height': height}
        if source == 'retina':
            self._retina_heatmap(**keyArgs)
        else:
            self._matr_heatmap(**keyArgs)
    
    def _retina_heatmap(self, normalize=1, title='', dist='bray-curtis', clust='ward', width=700, height=600):
        matrix = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
        # run our own R code
        matrix_file = Ipy.TMP_DIR+'/matrix.'+random_str()+'.tab'
        col_file = Ipy.TMP_DIR+'/col_clust.'+random_str()+'.txt'
        row_file = Ipy.TMP_DIR+'/row_clust.'+random_str()+'.txt'
        self.dump(fname=matrix_file, fformat='tab', normalize=normalize)
        rcmd = 'source("%s")\nMGRAST_dendrograms(file_in="%s", file_out_column="%s", file_out_row="%s", dist_method="%s", clust_method="%s", produce_figures="FALSE")\n'%(Ipy.LIB_DIR+'/dendrogram.r', matrix_file, col_file, row_file, dist, clust)
        ro.r(rcmd)
        cord, cdist = ordered_distance_from_file(col_file)
        rord, rdist = ordered_distance_from_file(row_file)
        data = { 'columns': self.ids(),
                 'rows': self.annotations(),
                 'colindex': cord,
                 'rowindex': rord,
                 'coldend': cdist,
                 'rowdend': rdist,
                 'data': matrix }
        keyArgs = { 'data': data, 'width': width, 'height': height }
        try:
            Ipy.RETINA.heatmap(**keyArgs)
        except:
            sys.stderr.write("Error producing heatmap\n")
    
    def _matr_heatmap(self, normalize=1, title='', dist='bray-curtis', clust='ward', width=700, height=600):
        matrix = self.NRmatrix if normalize and self.NRmatrix else self.Rmatrix
        fname  = Ipy.IMG_DIR+'/heatmap_'+random_str()+'.svg'
        labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        if not matrix:
            return None
        keyArgs = { 'labCol': ro.StrVector(labels),
                    'labRow': '',
                    'main': title,
                    'cexCol': 0.95,
                    'margins': ro.r.c(8,1) }
        ro.r.svg(fname)
        ro.r.heatmap(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def plot_annotation(self, normalize=1, ptype='row', width=800, height=800, x_rotate='0', title="", legend=True, subset=None, submg=None, onclick=None):
        matrix = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
        if not matrix:
            sys.stderr.write("Error producing chart: empty matrix\n")
            return
        if (not submg) or (len(submg) == 0):
            # default is all
            submg = self.ids()
        colors = google_palette(len(submg))
        labels = []
        data   = []
        # set retina data
        for i, col in enumerate(self.biom['columns']):
            if col['id'] in submg:
                data.append({'name': col['id'], 'data': [], 'fill': colors[i]})
        # populate data from matrix
        for r, row in enumerate(matrix):
            # only use subset rows if subset given
            if (subset is not None) and (self.biom['rows'][r]['id'] not in subset):
                continue
            labels.append(self.biom['rows'][r]['id'])
            for c, col in enumerate(self.biom['columns']):
                # only use submg cols
                if col['id'] in submg:
                    data[c]['data'].append(row[c])
        
        keyArgs = { 'btype': ptype,
                    'width': width,
                    'height': height,
                    'x_labels': labels,
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
            sys.stderr.write("Error producing chart\n")
    
    def _normalize_matrix(self):
        try:
            # can matr do it ?
            self.NRmatrix = ro.r.normalize(self.Rmatrix)
            self.NDmatrix = rMatrix_to_pyMatrix(self.NRmatrix, self.numAnnotations, self.numIDs)
        except:
            try:
                # run our own R code
                raw_file  = Ipy.TMP_DIR+'/raw.'+random_str()+'.tab'
                norm_file = Ipy.TMP_DIR+'/norm.'+random_str()+'.tab'
                self.dump(fname=raw_file, fformat='tab', normalize=0)
                rcmd = 'source("%s")\nMGRAST_preprocessing(file_in="%s", file_out="%s", produce_fig="FALSE")\n'%(Ipy.LIB_DIR+'/preprocessing.r', raw_file, norm_file)
                ro.r(rcmd)
                self.NDmatrix = matrix_from_file(norm_file, has_col_names=True, has_row_names=True)
                self.NRmatrix = pyMatrix_to_rMatrix(self.NDmatrix, self.numAnnotations, self.numIDs, normalize=1)
            except:
                sys.stderr.write("Error normalizing matrix (%s)\n"%self.id)

    def _dense_matrix(self):
        if not self.biom:
            return []
        if self.biom['matrix_type'] == 'dense':
            return self.matrix
        else:
            return sparse_to_dense(self.matrix, self.numAnnotations, self.numIDs)
