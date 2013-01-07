#!/usr/bin/env python

import sys, traceback
import retina
from ipyTools import *

class Metagenome:
    def __init__(self, mgid, metadata=True, stats=True, auth=None, def_name=None):
        self._retina = retina.Retina()
        self._auth = auth
        metagenome = self._get_metagenome(mgid, metadata)
        if metagenome is not None:
            for key, val in metagenome.iteritems():
                setattr(self, key, val)
        else:
            self.id = mgid
            self.name = None
        # get stats
        self.stats = None
        if stats:
            self._set_statistics()
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name
        
    def _get_metagenome(self, mgid, metadata):
        verb = 'full' if metadata else 'verbose'
        auth = '&auth='+self._auth if self._auth else ''
        return obj_from_url(API_URL+'metagenome/'+mgid+'?verbosity='+verb+auth)

    def _set_statistics(self):
        auth = '&auth='+self._auth if self._auth else ''
        self.stats = obj_from_url(API_URL+'metagenome_statistics/'+self.id+'?verbosity=full'+auth)
    
    def plot_taxon(self, ptype='pie', level='domain', parent=None):
        children = get_taxonomy(level, parent) if parent is not None else None
        self._plot_annotation('taxonomy', ptype, level, children)
    
    def plot_function(self, ptype='pie', source='Subsystems'):
        self._plot_annotation('ontology', ptype, source)
    
    def _plot_annotation(self, atype, ptype, level, names=None):
        if self.stats is None:
            self._set_statistics()
        data = []
        try:
            colors = google_palette(len(self.stats[atype][level]))
            for i, d in enumerate(self.stats[atype][level]):
                if (names is not None) and (d[0] not in names):
                    continue
                data.append({'name': d[0], 'data': [int(d[1])], 'fill': colors[i]})

            keyArgs = { 'btype': ptype,
                        'width': 700,
                        'height': 350,
                        'x_labels': '[""]',
                        'title': self.id+" "+level,
                        'target': self.id+"_"+level+'_'+random_str(),
                        'show_legend': True,
                        'legend_position': 'right',
                        'data': data }
            if atype == 'taxonomy':
                keyArgs['onclick'] = "'%s.plot_taxon(level=\"%s\", parent=\"'+params['series']+'\")'"%(self.defined_name, child_tax_level(level))
            self._retina.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            return None
