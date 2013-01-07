#!/usr/bin/env python

import sys, traceback
import analysis, retina
from ipyTools import *

class Project:
    def __init__(self, pid, metadata=True, auth=None, def_name=None):
        self._retina = retina.Retina()
        self._auth = auth
        project = self._get_project(pid, metadata)
        if project is not None:
            for key, val in project.iteritems():
                setattr(self, key, val)
        else:
            self.id = pid
            self.name = None
        self.stats = None
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name
    
    def _get_project(self, pid, metadata):
        verb = 'full' if metadata else 'verbose'
        auth = '&auth='+self._auth if self._auth else ''
        return obj_from_url(API_URL+'project/'+pid+'?verbosity='+verb+auth)
    
    def _set_statistics(self):
        self.stats = {}
        auth = '&auth='+self._auth if self._auth else ''
        for mgid in self.metagenomes():
            self.stats[mgid] = obj_from_url(API_URL+'metagenome_statistics/'+mgid+'?verbosity=verbose'+auth)
    
    def metagenomes(self):
        mlist = []
        if hasattr(self, 'analyzed'):
            mlist = map(lambda x: x[0], self.analyzed)
        return mlist
    
    def analysis_matrix(self, annotation='organism', level=None, resultType=None, source=None):
        return analysis.Analysis([self.id], annotation, level, resultType, source)

    def plot_taxon(self, ptype='column', level='domain', parent=None):
        children = get_taxonomy(level, parent) if parent is not None else None
        self._plot_annotation('taxonomy', ptype, level, children)

    def plot_function(self, ptype='column', source='Subsystems'):
        self._plot_annotation('ontology', ptype, source)

    def _plot_annotation(self, atype, ptype, level, names=None):
        if self.stats is None:
            self._set_statistics()
        data = []
        annD = {}
        try:
            colors = google_palette(len(self.stats))
            for i, mg, stats in enumerate(self.stats.iteritems()):
                for d in stats[atype][level]:
                    if (names is not None) and (d[0] not in names):
                        continue
                    data.append({'name': mg, 'data': [], 'fill': colors[i]})
                    annD[ d[0] ] = 1
            annL = annD.keys().sort()
            for d in data:
                annMG = {}
                for a, v in self.stats[d['name']][atype][level]:
                    annMG[a] = v
                for a in annL:
                    if a in annMG:
                        d['data'].append(int(annMG[a]))
                    else:
                        d['data'].append(0)
            
            keyArgs = { 'btype': ptype,
                        'width': 700,
                        'height': 350,
                        'x_labels': annL,
                        'title': self.id+" "+level,
                        'target': self.id+"_"+level+'_'+random_str(),
                        'show_legend': True,
                        'legend_position': 'right',
                        'data': data }
            print keyArgs
            self._retina.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            return None

