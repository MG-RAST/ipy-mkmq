#!/usr/bin/env python

import sys, os, json, traceback
from ipyTools import *

class Metagenome(object):
    def __init__(self, mgid, metadata=True, stats=True, auth=None, def_name=None, cache=None, mfile=None):
        self._auth = auth
        metagenome = None
        if mfile and os.path.isfile(mfile):
            # try load from file if given
            try:
                metagenome = json.load(open(mfile, 'rU'))
                sys.stdout.write("metagenome %s loaded from cache (%s)\n"%(mgid, cache))
            except:
                pass
        if metagenome is None:
            # load from api
            metagenome = self._get_metagenome(mgid, metadata)
            if cache and metagenome and os.path.isdir(cache):
                # cache it if dir given and not loaded from file
                try:
                    json.dump(metagenome, open(cache+'/'+mgid+'.json', 'rU'))
                    sys.stdout.write("metagenome %s saved to cache (%s)\n"%(mgid, cache))
                except:
                    pass
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
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                def_name = None
        self.defined_name = def_name
        
    def _get_metagenome(self, mgid, metadata):
        verb = 'full' if metadata else 'verbose'
        auth = '&auth='+self._auth if self._auth else ''
        return obj_from_url(Ipy.API_URL+'metagenome/'+mgid+'?verbosity='+verb+auth)

    def _set_statistics(self):
        auth = '&auth='+self._auth if self._auth else ''
        self.stats = obj_from_url(Ipy.API_URL+'metagenome_statistics/'+self.id+'?verbosity=full'+auth)
    
    def show_metadata(self):
        mdTable = []
        if hasattr(self, 'metadata'):
            for cat, data in self.metadata.iteritems():
                for field, value in data['data'].iteritems():
                    mdTable.append([cat, field, value])
        if len(mdTable) == 0:
            sys.stderr.write("No metadata to display\n")
        keyArgs = { 'width': 700,
                    'height': 600,
                    'target': self.id+"_metadata_"+random_str(),
                    'data': {'data': mdTable, 'header': ['category', 'field', 'value']},
                    'rows_per_page': 20 }
        if Ipy.DEBUG:
            print keyArgs
        try:
            Ipy.RETINA.table(**keyArgs)
        except:
            sys.stderr.write("Error producing metadata table\n")
    
    def piechart_taxon(self, level='domain', parent=None, title=''):
        children = get_taxonomy(level, parent) if parent is not None else None
        self._piechart('taxonomy', level, names=children, title=title)
    
    def piechart_function(self, source='Subsystems', title=''):
        self._piechart('ontology', source, title=title)
    
    def _piechart(self, atype, level, names=None, title=''):
        if self.stats is None:
            self._set_statistics()
        data = []
        try:
            colors = google_palette(len(self.stats[atype][level]))
            for i, d in enumerate(self.stats[atype][level]):
                if (names is not None) and (d[0] not in names):
                    continue
                data.append({'name': d[0], 'data': [int(d[1])], 'fill': colors[i]})
            lheight = len(self.stats[atype][level])*30
            lwidth  = int(len(max(self.stats[atype][level], key=len))*7.2)
            keyArgs = { 'btype': 'pie',
                        'width': 700 + int((float(lwidth)/2)),
                        'height': 350,
                        'x_labels': [""],
                        'title': title,
                        'target': self.id+"_"+level+'_'+random_str(),
                        'show_legend': True,
                        'legendArea': [0.80, 0.05, lwidth, lheight],
                        'data': data }
            if atype == 'taxonomy':
                keyArgs['onclick'] = "'%s.plot_taxon(level=\"%s\", parent=\"'+params['series']+'\")'"%(self.defined_name, child_level(level, htype='taxonomy'))
            if Ipy.DEBUG:
                print keyArgs
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing %s chart"%atype)
