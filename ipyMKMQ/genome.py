#!/usr/bin/env python

import sys, os, json, traceback
from ipyTools import *

class Genome(object):
    """Class representation of Genome object:
        "id"       : [ 'string',  'unique object identifier' ],
        "name"     : [ 'string',  'human readable identifier' ],
        "metadata" : [ 'hash',    'key value pairs describing metadata' ],
        "version"  : [ 'integer', 'version of the object' ],
        "url"      : [ 'uri',     'resource location of this object instance' ],
        "contigs"  : [ 'hash', {
                        "id" :     [ 'string',  'contig ID' ],
                        "length" : [ 'integer', 'contig length' ],
                       }]
    """
    def __init__(self, gid, metadata=True, stats=True, auth=None, def_name=None, cache=None, mfile=None):
        self._auth = auth
        genome = None
        if mfile and os.path.isfile(mfile):
            # try load from file if given
            try:
                genome = json.load(open(mfile, 'rU'))
                if Ipy.DEBUG:
                    sys.stdout.write("genome %s loaded from cache (%s)\n"%(gid, cache))
            except:
                pass
        if genome is None:
            # load from api
            genome = self._get_genome(gid, metadata)
            if cache and genome and os.path.isdir(cache):
                # cache it if dir given and not loaded from file
                try:
                    json.dump(genome, open(cache+'/'+gid+'.json', 'w'))
                    if Ipy.DEBUG:
                        sys.stdout.write("genome %s saved to cache (%s)\n"%(gid, cache))
                except:
                    pass
        if genome is not None:
            for key, val in genome.iteritems():
                setattr(self, key, val)
        else:
            self.id = gid
            self.name = None
        # hack to get variable name
        if def_name == None:
            try:
                (filename, line_number, function_name, text) = traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self.defined_name = def_name
        
    def _get_genome(self, gid, metadata):
        verb = 'full' if metadata else 'verbose'
        return obj_from_url(Ipy.API_URL+'genome/'+gid, self._auth)
    
    def show_metadata(self):
        mdTable = []
        if hasattr(self, 'metadata'):
            for cat, data in self.metadata.iteritems():
                for field, value in data['data'].iteritems():
                    mdTable.append([cat, field, value])
        if len(mdTable) == 0:
            sys.stderr.write("No metadata to display\n")
        keyArgs = {
            'width': 700,
            'height': 600,
            'target': self.id+"_metadata_"+random_str(),
            'data': {
                'data': mdTable,
                'header': ['category', 'field', 'value']
            },
            'rows_per_page': 20
        }
        if Ipy.DEBUG:
            print keyArgs
        try:
            Ipy.RETINA.table(**keyArgs)
        except:
            sys.stderr.write("Error producing metadata table\n")
    
    def barchart_contigs(self, title=''):
        contigs = get_contigs()
        self._barchart('Contigs', contigs, title=title)
    
    def _barchart(self, data, names=None, title=''):
        try:
            colors = google_palette(len(data))
            for i, d in enumerate(self.stats[atype][level]):
                if (names is not None) and (d[0] not in names):
                    continue
                data.append({'name': d[0], 'data': [int(d[1])], 'fill': colors[i]})
            lheight = len(self.stats[atype][level])*30
            lwidth  = int(len(max(self.stats[atype][level], key=len))*7.2)
            keyArgs = { 'btype': 'bar',
                        'width': 700 + int((float(lwidth)/2)),
                        'height': 350,
                        'x_labels': [""],
                        'title': title,
                        'target': self.id+"_"+level+'_'+random_str(),
                        'show_legend': True,
                        'legendArea': [0.80, 0.05, lwidth, lheight],
                        'data': data }
            if Ipy.DEBUG:
                print keyArgs
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing %s chart"%atype)
