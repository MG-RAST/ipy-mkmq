#!/usr/bin/env python

import shutil, hashlib
import math, urllib, sys, os, traceback
import rpy2.robjects as ro
from metagenome import Metagenome
from ipyTools import *
from collections import defaultdict
from datetime import datetime

class AnalysisSet(object):
    """Class for working with a set of Analysis objects:
        - Creates an Analysis object for each taxonimic level and functional level
        - allows barchart and heatmap navigation through hierarchies (drilldowns)
        - caches data locally for fast re-analysis
    """
    def __init__(self, ids=[], auth=None, method='WGS', all_values=False, cache=None, def_name=None, reset_cache=False):
        if cache is None:
            cache = random_str()
        self.method = method
        self.cache = cache
        self._path = Ipy.NB_DIR+'/'+cache
        self._auth = auth
        self.all_mgs = ids
        self.display_mgs = self.all_mgs
        tax_source = ''
        if self.method == 'WGS':
            tax_source = 'M5NR'
        elif self.method == 'Amplicon':
            tax_source = 'M5RNA'
        else:
            sys.stderr.write("Error: invalid method (%s), use one of 'WGS' or 'Amplicon'"%self.method)
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name
        # set cache
        self.cache_time = self._set_cache(reset_cache)
        # get data
        for tax in Ipy.TAX_SET:
            values = {}
            if all_values:
                for val in Ipy.VALUES:
                    values[val] = self._get_analysis(ids, 'organism', tax, val, tax_source)
            else:
                values['abundance'] = self._get_analysis(ids, 'organism', tax, 'abundance', tax_source)
            setattr(self, tax, values)
                
        if self.method == 'WGS':
            for ont in Ipy.ONT_SET:
                values = {}
                if all_values:
                    for val in Ipy.VALUES:
                        values[val] = self._get_analysis(ids, 'function', ont, val, 'Subsystems')
                else:
                    values['abundance'] = self._get_analysis(ids, 'function', ont, 'abundance', 'Subsystems')
                setattr(self, ont, values)

    def set_display_mgs(self, ids=[]):
        if (not ids) or (len(ids) == 0):
            sys.stdout.write("setting %s.display_mgs to all metagenomes in set\n"%self.defined_name)
            self.display_mgs = self.all_mgs
        else:
            self.display_mgs = ids

    def _set_cache(self, reset=False):
        # force re-cache
        if reset and os.path.isdir(self._path):
            shutil.rmtree(self._path)
        # test if exists
        if os.path.isdir(self._path) and os.path.isfile(self._path+'/CACHE_TIME'):
            thdl  = open(self._path+'/CACHE_TIME', 'rU')
            ctime = thdl.read().strip()
            thdl.close()
            sys.stdout.write("analysis-set '%s' loaded from cache %s (%s)\n"%(self.defined_name, self.cache, ctime))
            return ctime
        # set dir and time
        os.mkdir(self._path)
        ctime = self._timestamp_cache()
        sys.stdout.write("analysis-set '%s' saved to cache %s (%s)\n"%(self.defined_name, self.cache, ctime))
        return ctime

    def _timestamp_cache(self, ctime=None):
        if os.path.isdir(self._path):
            ctime = ctime if ctime else str(datetime.now())
            thdl  = open(self._path+'/CACHE_TIME', 'w')
            thdl.write(ctime+"\n")
            thdl.close()
            return ctime
        else:
            return None

    def dump(self, force=False):
        if not os.path.isdir(self._path):
            self.cache_time = self._set_cache()
        # dump individual files
        for key in Ipy.TAX_SET:
            item = getattr(self, key)
            for analysis in item.itervalues():
                fname = self._path+'/'+hashlib.md5(analysis.id).hexdigest()+'.biom'
                if force or (not os.path.isfile(fname)):
                    analysis.dump(fname=fname, fformat='biom')
        if self.method == 'WGS':
            for key in Ipy.ONT_SET:
                item = getattr(self, key)
                for analysis in item.itervalues():
                    fname = self._path+'/'+hashlib.md5(analysis.id).hexdigest()+'.biom'
                    if force or (not os.path.isfile(fname)):
                        analysis.dump(fname=fname, fformat='biom')
        if force:
            self.cache_time = self._timestamp_cache()

    def _get_analysis(self, ids, annotation, level, result_type, source):
        # this needs to be created same way as matrix api builds it
        matrix_id = "_".join(sorted(ids))+"_"+"_".join([annotation, level, source, result_type])
        matrix_id += "_%d_%d_%d"%(Ipy.MATRIX['e_val'], Ipy.MATRIX['ident'], Ipy.MATRIX['alen'])
        # load from client cache if exists
        matrix_md5 = hashlib.md5(matrix_id).hexdigest()
        biom_file  = self._path+'/'+matrix_md5+'.biom'
        if os.path.isfile(biom_file):
            if Ipy.DEBUG:
                sys.stdout.write("loading %s.biom (%s) from cache %s ... \n"%(matrix_md5, matrix_id, self.cache))
            return Analysis(bfile=biom_file, auth=self._auth)
        else:
            if Ipy.DEBUG:
                sys.stdout.write("loading %s through api ... \n"%matrix_id)
            keyArgs = Ipy.MATRIX
            keyArgs['ids'] = ids
            keyArgs['annotation'] = annotation
            keyArgs['level'] = level
            keyArgs['result_type'] = result_type
            keyArgs['source'] = source
            if self._auth:
                keyArgs['auth'] = self._auth
            thisAnalysis = Analysis(**keyArgs)
            thisAnalysis.dump(fname=biom_file, fformat='biom')
            if Ipy.DEBUG:
                sys.stdout.write("caching %s.biom (%s) to %s ... \n"%(matrix_md5, matrix_id, self.cache))
            return thisAnalysis
    
    def barchart(self, annot='organism', level='domain', parent=None, width=800, height=0, title="", legend=True, normalize=1):
        children = []
        if parent and (len(parent) > 0):
            for p in parent:
                children.extend( get_hierarchy(htype=annot, level=level, parent=p) )
        if children and (len(children) > 0):
            children = filter(lambda x: x, children)
        keyArgs = { 'normalize': normalize,
                    'width': width,
                    'height': height,
                    'x_rotate': '0',
                    'title': title,
                    'legend': legend,
                    'submg': self.display_mgs,
                    'subset': children }
        if child_level(level, htype=annot):
            click_opts = (self.defined_name, child_level(level, htype=annot), annot, normalize, width, height, title, 'True' if legend else 'False')
            keyArgs['onclick'] = "'%s.barchart(level=\"%s\", parent=[\"'+params['label']+'\"], annot=\"%s\", normalize=%d, width=%d, height=%d, title=\"%s\", legend=%s)'"%click_opts
        if Ipy.DEBUG:
            print annot, level, child_level(level, htype=annot), keyArgs
        to_plot = getattr(self, level)
        to_plot['abundance'].barchart(**keyArgs)
        
    def heatmap(self, annot='organism', level='domain', parent=None, width=700, height=600, normalize=1, dist='bray-curtis', clust='ward'):
        children = []
        if parent and (len(parent) > 0):
            for p in parent:
                children.extend( get_hierarchy(htype=annot, level=level, parent=p) )
        if children and (len(children) > 0):
            children = filter(lambda x: x, children)
        keyArgs = { 'normalize': normalize,
                    'width': width,
                    'height': height,
                    'dist': dist,
                    'clust': clust,
                    'submg': self.display_mgs,
                    'subset': children }
        if child_level(level, htype=annot):
            click_opts = (self.defined_name, child_level(level, htype=annot), annot, normalize, width, height, dist, clust)
            keyArgs['onclick'] = "'%s.heatmap(level=\"%s\", parent=\"'+sel_names+'\", annot=\"%s\", normalize=%d, width=%d, height=%d, dist=\"%s\", clust=\"%s\")'"%click_opts
        if Ipy.DEBUG:
            print annot, level, keyArgs
        to_plot = getattr(self, level)
        to_plot['abundance']._retina_heatmap(**keyArgs)
        

class Analysis(object):
    """Class representation of Matrix object:
        self.biom (BIOM format):
            "id"                   : [ 'string', 'unique object identifier' ],
            "format"               : [ 'string', 'format specification name' ],
            "format_url"           : [ 'string', 'url to the format specification' ],
            "type"                 : [ 'string', 'type of the data in the return table (taxon, function or gene)' ],
            "generated_by"         : [ 'string', 'identifier of the data generator' ],
            "date"                 : [ 'date',   'time the output data was generated' ],
            "matrix_type"          : [ 'string', 'type of the data encoding matrix (dense or sparse)' ],
            "matrix_element_type"  : [ 'string', 'data type of the elements in the return matrix' ],
            "matrix_element_value" : [ 'string', 'result_type of the elements in the return matrix' ],
            "shape"                : [ 'list', ['integer', 'list of the dimension sizes of the return matrix'] ],
            "rows"                 : [ 'list', ['object', [{'id'       => ['string', 'unique annotation text'],
                                                            'metadata' => ['hash', 'key value pairs describing metadata']}, "rows object"]] ],
            "columns"              : [ 'list', ['object', [{'id'       => ['string', 'unique metagenome identifier'],
            	                                            'metadata' => ['hash', 'key value pairs describing metadata']}, "columns object"]] ],
            "data"                 : [ 'list', ['list', ['float', 'the matrix values']] ]
        self.id       : BIOM id
        self.numIDs   : BIOM column count
        self.numAnnot : BIOM row count
        self.Dmatrix  : dense matrix of BIOM data
        self.Rmatrix  : R-format dense matrix
        self.NDmatrix : normalized dense matrix
        self.NRmatrix : normalized R-format dense matrix
        
        Visualizations:
            self.dump()     : produce file or string of BIOM or tab-deliminated matrix
            self.boxplot()  : boxplot display
            self.barchart() : horizontal barchart of metagenomes / annotations
            self.pco()      : pco plot of metagenomes
            self.heatmap()  : dendogram of metagenomes / annotations
    """
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
        self.numIDs = self.biom['shape'][1] if self.biom else 0
        self.numAnnot = self.biom['shape'][0] if self.biom else 0
        self.Dmatrix  = self._dense_matrix()  # count dense matrix
        self.Rmatrix  = pyMatrix_to_rMatrix(self.Dmatrix, self.numAnnot, self.numIDs) # R count matrix object
        self.NDmatrix = None  # normalized dense matrix
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

    def find_annotation(self, text):
        if not self.biom:
            return []
        str_re = re.compile(text, re.IGNORECASE)
        annot = set()
        hier = ''
        if self.biom['type'].startswith('Taxon'):
            hier = 'taxonomy'
        elif self.biom['type'].startswith('Function'):
            hier = 'ontology'
        for r in self.biom['rows']:
            if str_re.search(r['id']):
                annot.add(r['id'])
            elif r['metadata'] and hier and (hier in r['metadata']) and str_re.search(r['metadata'][hier][-1]):
                annot.add(r['metadata'][hier][-1])
        return list(annot)

    def sub_matrix(self, normalize=0, cols=None, rows=None):
        all_annot = self.annotations()
        all_mgids = self.ids()
        matrix = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
        # validate rows
        aIndex = []
        if rows and (len(rows) > 0):
            aIndex = range(len(all_annot))
        else:
            for r in rows:
                try:
                    a = all_annot.index(r)
                except (ValueError, AttributeError):
                    continue
                aIndex.append(a)
        # validate cols
        mgids  = []
        mIndex = []
        if cols and (len(cols) > 0):
            mgids  = all_mgids
            mIndex = range(len(all_mgids))
        else:
            for c in cols:
                try:
                    m = all_mgids.index(c)
                except (ValueError, AttributeError):
                    continue
                mgids.append(all_mgids[m])
                mIndex.append(m)
        # build matrix
        sub_matrix = []
        good_rows = []
        for i in aIndex:
            rdata = []
            for j in mIndex:
                rdata.append(matrix[i][j])
            if sum(rdata) > 0:
                sub_matrix.append(rdata)
                good_rows.append(all_annot[i])
        return good_rows, mgids, matrix

    def dump(self, fname=None, fformat='biom', normalize=0, rows=None, cols=None, metadata=False):
        output = ""
        if not self.biom:
            sys.stderr.write("Error dumping %s, no data\n"%self.id)
            return
        if fformat == 'biom':
            # biom dump
            output = json.dumps(self.biom)
        else:
            # tab deleminted dump / option for sub-dump
            all_annot = self.annotations(metadata=metadata)
            all_mgids = self.ids()
            matrix = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
            annot  = rows if rows and (len(rows) > 0) else all_annot
            mgids  = cols if cols and (len(cols) > 0) else all_mgids
            output = "\t%s\n"%"\t".join(mgids)
            for a in annot:
                try:
                    r = all_annot.index(a)
                except:
                    sys.stderr.write("Error: '%s' is not in annotations of %s"%(a, self.id))
                    return None
                output += a
                for m in mgids:
                    try:
                        c = all_mgids.index(m)
                    except:
                        sys.stderr.write("Error: '%s' is not in metagenomes of %s"%(m, self.id))
                        return None
                    output += "\t"+str(matrix[r][c])
                output += "\n"
        if fname:
            open(fname, 'w').write(output)
        else:
            return output

    def ids(self):
        if not self.biom:
            return None
        return map(lambda x: x['id'], self.biom['columns'])
    
    def names(self):
        if not self.biom:
            return None
        return map(lambda x: x['name'], self.biom['columns'])
    
    def annotations(self, metadata=False):
        if not self.biom:
            return None
        if not metadata:
            return map(lambda x: x['id'], self.biom['rows'])
        else:
            ann = []
            for r in self.biom['rows']:
                if ('metadata' in r) and ('ontology' in r['metadata']):
                    ann.append(r['metadata']['ontology'][-1])
                else:
                    ann.append(r['id'])
            return ann
    
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
    
    def boxplot(self, normalize=1, title=''):
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
        if Ipy.DEBUG:
            print fname, keyArgs
        ro.r.svg(fname)
        ro.r.boxplot(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname
    
    def pco(self, normalize=1, dist='bray-curtis', title=''):
        matrix = self.NRmatrix if normalize and self.NRmatrix else self.Rmatrix
        fname  = Ipy.IMG_DIR+'/pco_'+random_str()+'.svg'
        labels = map(lambda x: x['id']+"\n"+x['name'], self.biom['columns'])
        if not matrix:
            return None
        keyArgs = { 'labels': ro.StrVector(labels),
                    'main': title,
                    'method': dist,
                    'comp': ro.r.c(1,2,3) }
        if Ipy.DEBUG:
            print fname, keyArgs
        ro.r.svg(fname)
        ro.r.pco(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname

    def heatmap(self, normalize=1, title='', dist='bray-curtis', clust='ward', width=700, height=600, source='retina'):
        if source == 'retina':
            self._retina_heatmap(normalize=normalize, dist=dist, clust=clust, width=width, height=height)
            return
        else:
            return self._matr_heatmap(normalize=normalize, title=title)
    
    def _retina_heatmap(self, normalize=1, dist='bray-curtis', clust='ward', width=700, height=600, submg=None, subset=None):
        # default is all
        all_annot = self.annotations(metadata=True)
        if (not submg) or (len(submg) == 0):
            submg = self.ids()
        if (not subset) or (len(subset) == 0):
            subset = all_annot
        else:
            subset = filter(lambda x: x in all_annot, subset)
        # run our own R code
        matrix_file = Ipy.TMP_DIR+'/matrix.'+random_str()+'.tab'
        col_file = Ipy.TMP_DIR+'/col_clust.'+random_str()+'.txt'
        row_file = Ipy.TMP_DIR+'/row_clust.'+random_str()+'.txt'
        self.dump(fname=matrix_file, fformat='tab', normalize=normalize, rows=subset, cols=submg, metadata=True)
        rcmd = 'source("%s")\nMGRAST_dendrograms(file_in="%s", file_out_column="%s", file_out_row="%s", dist_method="%s", clust_method="%s", produce_figures="FALSE")\n'%(Ipy.LIB_DIR+'/dendrogram.r', matrix_file, col_file, row_file, dist, clust)
        ro.r(rcmd)
        cord, cdist = ordered_distance_from_file(col_file)
        rord, rdist = ordered_distance_from_file(row_file)
        sub_matrix  = matrix_from_file(matrix_file)
        data = { 'columns': submg,
                 'rows': subset,
                 'colindex': cord,
                 'rowindex': rord,
                 'coldend': cdist,
                 'rowdend': rdist,
                 'data': sub_matrix }
        lwidth  = len(max(subset, key=len)) * 7.2
        keyArgs = { 'data': data,
                    'width': int(width+lwidth),
                    'height': height,
                    'target': 'div_heatmap_'+random_str(),
                    'tree_width': 200,
                    'legend_width': int(lwidth) }
        if Ipy.DEBUG:
            print keyArgs
        try:
            Ipy.RETINA.heatmap(**keyArgs)
        except:
            sys.stderr.write("Error producing heatmap\n")

    def _matr_heatmap(self, normalize=1, title=''):
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
        if Ipy.DEBUG:
            print fname
        ro.r.svg(fname)
        ro.r.heatmap(matrix, **keyArgs)
        ro.r("dev.off()")
        return fname

    def barchart(self, normalize=1, width=800, height=0, x_rotate='0', title="", legend=True, subset=None, submg=None, onclick=None):
        matrix  = self.NDmatrix if normalize and self.NDmatrix else self.Dmatrix
        all_ids = self.ids()
        if not matrix:
            sys.stderr.write("Error producing chart: empty matrix\n")
            return
        if (not submg) or (len(submg) == 0):
            # default is all
            submg = all_ids
        colors = google_palette(len(submg))
        labels = []
        data   = []
        # set retina data
        for i, m in enumerate(submg):
            if m in all_ids:
                data.append({'name': m, 'data': [], 'fill': colors[i]})
        # populate data from matrix
        for r, row in enumerate(matrix):
            # only use subset rows if subset given
            rdata = self.biom['rows'][r]
            if subset and (len(subset) > 0):
                if (rdata['id'] not in subset) or (('ontology' in rdata['metadata']) and (rdata['metadata']['ontology'][-1] not in subset)):
                    continue
            rname = rdata['metadata']['ontology'][-1] if 'ontology' in rdata['metadata'] else rdata['id']
            labels.append(rname)
            # only use submg cols
            for i, m in enumerate(submg):
                c = all_ids.index(m)
                data[i]['data'].append(toNum(row[c]))
        height  = height if height else len(labels)*len(submg)*7.5
        lheight = min(height, len(submg)*35)
        lwidth  = len(max(labels, key=len)) * 7.2
        cwidth  = 0.85 if legend else 0.99
        keyArgs = { 'btype': 'row',
                    'width': int(width+lwidth),
                    'height': int(height),
                    'x_labels': labels,
                    'x_labels_rotation': x_rotate,
                    'title': title,
                    'target': 'div_graph_'+random_str(),
                    'show_legend': legend,
                    'legendArea': [0.87, 0.05, 0.2, int(lheight)],
                    'chartArea': [int(lwidth), 0.02, cwidth, 0.95],
                    'data': data,
                    'onclick': onclick }
        if normalize and self.NDmatrix:
            keyArgs['y_labeled_tick_interval'] = 0.1
        if Ipy.DEBUG:
            print submg, subset, keyArgs
        try:
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart\n")

    def _normalize_matrix(self):
        try:
            # can matr do it ?
            self.NRmatrix = ro.r.normalize(self.Rmatrix)
            self.NDmatrix = rMatrix_to_pyMatrix(self.NRmatrix, self.numAnnot, self.numIDs)
        except:
            try:
                # run our own R code
                raw_file  = Ipy.TMP_DIR+'/raw.'+random_str()+'.tab'
                norm_file = Ipy.TMP_DIR+'/norm.'+random_str()+'.tab'
                self.dump(fname=raw_file, fformat='tab', normalize=0)
                rcmd = 'source("%s")\nMGRAST_preprocessing(file_in="%s", file_out="%s", produce_fig="FALSE")\n'%(Ipy.LIB_DIR+'/preprocessing.r', raw_file, norm_file)
                ro.r(rcmd)
                self.NDmatrix = matrix_from_file(norm_file, has_col_names=True, has_row_names=True)
                self.NRmatrix = pyMatrix_to_rMatrix(self.NDmatrix, self.numAnnot, self.numIDs, normalize=1)
            except:
                sys.stderr.write("Error normalizing matrix (%s)\n"%self.id)

    def _dense_matrix(self):
        if not self.biom:
            return []
        if self.biom['matrix_type'] == 'dense':
            return self.biom['data']
        else:
            return sparse_to_dense(self.biom['data'], self.numAnnot, self.numIDs)
