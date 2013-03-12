#!/usr/bin/env python

import sys, os, json, traceback
from ipyTools import *

class Metagenome(object):
    """Class representation of Metagenome object:
        "id"       : [ 'string', 'unique object identifier' ],
        "name"     : [ 'string', 'human readable identifier' ],
        "library"  : [ 'reference library', 'reference to the related library object' ],
        "sample"   : [ 'reference sample',  'reference to the related sample object' ],
        "project"  : [ 'reference project', 'reference to the project object' ],
        "metadata" : [ 'hash',    'key value pairs describing metadata' ],
        "created"  : [ 'date',    'time the object was first created' ],
        "version"  : [ 'integer', 'version of the object' ],
        "url"      : [ 'uri',     'resource location of this object instance' ],
        "status"   : [ 'cv', [ ['public', 'object is public'],
        					   ['private', 'object is private'] ] ],
        "sequence_type" : [ 'string', 'sequencing type' ]
        "stats"      : "id" : [ 'string', 'unique metagenome id' ],
                       "length_histogram" : { "upload" : [ 'list', 'length distribution of uploaded sequences' ],
                                              "post_qc" : [ 'list', 'length distribution of post-qc sequences' ] },
                       "gc_histogram" : { "upload" : [ 'list', 'gc % distribution of uploaded sequences' ],
                                          "post_qc" : [ 'list', 'gc % distribution of post-qc sequences' ] },
                       "qc" : { "kmer" : { "6_mer"  : {"columns" : ['list', 'names of columns'], "data" : ['list', 'kmer 6 counts']},
                                           "15_mer" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'kmer 15 counts']} },
                                "drisee" : { "counts" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'drisee count profile']},
                                             "percents" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'drisee percent profile']},
                                             "summary" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'drisee summary stats']} },
                                "bp_profile" : { "counts" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'nucleotide count profile']},
                                                 "percents" : {"columns" : ['list', 'names of columns'], "data" : ['list', 'nucleotide percent profile']} }
                               },
                       "sequence_stats" : [ 'hash', 'statistics on sequence files of all pipeline stages' ],
                       "taxonomy" : { "species" : [ 'list', 'species counts' ],
                                      "genus" : [ 'list', 'genus counts' ],
                                      "family" : [ 'list', 'family counts' ],
                                      "order" : [ 'list', 'order counts' ],
                                      "class" : [ 'list', 'class counts' ],
                                      "phylum" : [ 'list', 'phylum counts' ],
                                      "domain" : [ 'list', 'domain counts' ] },
                       "ontology" : { "COG" : [ 'list', 'COG counts' ],
                                      "KO" : [ 'list', 'KO counts' ],
		                              "NOG" : [ 'list', 'NOG counts' ],
		                              "Subsystems" : [ 'list', 'Subsystem counts' ] },
                       "source" : [ 'hash', 'evalue and % identity counts per source' ],
	                   "rarefaction" : [ 'list', 'rarefaction coordinate data' ]
    """
    def __init__(self, mgid, stats=True, auth=None, def_name=None, cache=False):
        self._auth  = auth
        self._cfile = Ipy.CCH_DIR+'/'+mgid+'.json'
        self.stats  = None
        metagenome  = None
        if cache and os.path.isfile(self._cfile):
            # try load from cache if given
            try:
                metagenome = json.load(open(self._cfile, 'rU'))
                if Ipy.DEBUG:
                    sys.stdout.write("metagenome %s loaded from cached file (%s)\n"%(mgid, self._cfile))
            except:
                pass
        if metagenome is None:
            # load from api
            metagenome = self._get_metagenome(mgid)
            if metagenome and cache and os.path.isdir(Ipy.CCH_DIR):
                # save to cache if given
                try:
                    json.dump(metagenome, open(self._cfile, 'w'))
                    if Ipy.DEBUG:
                        sys.stdout.write("metagenome %s saved to cached file (%s)\n"%(mgid, self._cfile))
                except:
                    pass
        if metagenome is not None:
            for key, val in metagenome.iteritems():
                setattr(self, key, val)
        else:
            sys.stderr.write("ERROR: unable to load metagenome %s\n"%mgid)
            self.id = mgid
            self.name = None
            return
        # get stats
        if stats:
            self._set_statistics()
        # hack to get variable name
        if def_name == None:
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self._defined_name = def_name
        
    def _get_metagenome(self, mgid):
        if Ipy.DEBUG:
            sys.stdout.write("Loading metagenome %s from API ...\n"%mgid)
        return obj_from_url(Ipy.API_URL+'/metagenome/'+mgid+'?verbosity=full', self._auth)

    def _set_statistics(self):
        if Ipy.DEBUG:
            sys.stdout.write("Loading metagenome %s statistics from API ...\n"%self.id)
        self.stats = obj_from_url(Ipy.API_URL+'/metagenome_statistics/'+self.id+'?verbosity=full', self._auth)
    
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
                qname = self._defined_name.replace("'", "\\\'")
                keyArgs['onclick'] = '%s.piechart_taxon(level="%s", parent="\'+params[\'series\']+\'")'%(qname, child_level(level, htype='taxonomy'))
            if Ipy.DEBUG:
                print keyArgs
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing %s chart"%atype)
