#!/usr/bin/env python

import sys
import analysis
import retina
from ipyTools import *

RETINA = retina.Retina()

class Project:
    def __init__(self, pid, metadata=True):
        project = self._get_project(pid, metadata)
        if project is not None:
            for key, val in project.iteritems():
                self[key] = val
        else:
            self.id = pid
            self.name = None
        self.stats = None
    
    def _get_project(self, pid, metadata):
        verb = 'full' if metadata else 'verbose'
        return obj_from_url(API_URL+'project/'+pid+'?verbosity='+verb)
    
    def _set_statistics(self):
        self.stats = {}
        for mgid in self.metagenomes():
            self.stats[mgid] = obj_from_url(API_URL+'metagenome_statistics/'+mgid+'?verbosity=verbose')
    
    def metagenomes(self):
        mlist = []
        if 'analyzed' in self:
            mlist = map(lambda x: x[0], self.analyzed)
        return mlist
    
    def analysis_matrix(self, annotation='organism', level=None, resultType=None, source=None):
        return analysis.Analysis([self.id], annotation, level, resultType, source)

    def _plot_annotation(self, atype, ptype, level, names=None):
        if self.stats is None:
            self.stats = self._set_statistics()
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
                        d['data'].append(annMG[a])
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
                        'title_settings': {'font-size': '18px', 'font-weight': 'bold', 'x': 0, 'text-anchor': 'start'},
                        'data': data }
            RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            return None

