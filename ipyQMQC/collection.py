#!/usr/bin/env python

import json, sys, traceback
import analysis
from metagenome import Metagenome
from ipyTools import *

class Collection(object):
    def __init__(self, mgids, metadata=True, stats=True, auth=None, def_name=None):
        self._auth  = auth
        self._stats = stats
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name
        # get metagenomes
        self._mgids = mgids
        self.metagenomes = self._get_metagenomes(mgids, metadata, stats)
    
    def _get_metagenomes(self, mgids, metadata, stats):
        mgs = {}
        for mg in mgids:
            keyArgs = { 'metadata': metadata,
                        'stats': stats,
                        'auth': self._auth,
                        'def_name': '%s.metagenomes["%s"]'%(self.defined_name, mg)
                       }
            mgs[mg] = Metagenome(mg, **keyArgs)
        return mgs
    
    def _set_statistics(self):
        self._stats = True
        for mg in self.metagenomes.itervalues():
            mg._set_statistics()
    
    def metagenomes(self):
        return self._mgids
        
    def analysis_matrix(self, annotation='organism', level=None, resultType=None, source=None):
        keyArgs = { 'ids': self.metagenomes(),
                    'annotation': annotation,
                    'level': level,
                    'resultType': resultType,
                    'source': source,
                    'auth': self._auth }
        return analysis.Analysis(**keyArgs)

    def plot_taxon(self, ptype='row', level='domain', parent=None, width=900, height=400, x_rotate='0', title=None, legend=True):
        children = get_taxonomy(level, parent) if parent is not None else None
        self._plot_annotation('taxonomy', ptype, level, width, height, x_rotate, title, legend, names=children)

    def plot_function(self, ptype='row', source='Subsystems', width=900, height=400, x_rotate='0', title=None, legend=True):
        self._plot_annotation('ontology', ptype, source, width, height, x_rotate, title, legend)

    def _plot_annotation(self, atype, ptype, level, width, height, x_rotate, title, legend, names=None):
        if not self._stats:
            self._set_statistics()
        data = []
        annD = {}
        colors = google_palette(len(self.metagenomes))
        for i, item in enumerate(self.metagenomes.iteritems()):
            mid, mg = item
            data.append({'name': mid, 'data': [], 'fill': colors[i]})
            for d in mg.stats[atype][level]:
                if (names is not None) and (d[0] not in names):
                    continue
                annD[ d[0] ] = 1
        annL = sorted(annD.keys())
        for d in data:
            annMG = {}
            for a, v in self.metagenomes[d['name']].stats[atype][level]:
                annMG[a] = v
            for a in annL:
                if a in annMG:
                    d['data'].append(int(annMG[a]))
                else:
                    d['data'].append(0)
            
        keyArgs = { 'btype': ptype,
                    'width': width,
                    'height': height,
                    'x_labels': json.dumps(annL),
                    'x_labels_rotation': x_rotate,
                    'title': title,
                    'target': '_'.join(self.metagenomes())+"_"+level+'_'+random_str(),
                    'show_legend': legend,
                    'legend_position': 'right',
                    'data': data }
        try:
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            return None
    