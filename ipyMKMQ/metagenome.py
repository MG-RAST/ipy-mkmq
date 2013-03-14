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
        # set display
        self.display = MetagenomeDisplay(self, self._defined_name+'.display')
        
    def _get_metagenome(self, mgid):
        if Ipy.DEBUG:
            sys.stdout.write("Loading metagenome %s from API ...\n"%mgid)
        return obj_from_url(Ipy.API_URL+'/metagenome/'+mgid+'?verbosity=full', self._auth)

    def _set_statistics(self):
        if Ipy.DEBUG:
            sys.stdout.write("Loading metagenome %s statistics from API ...\n"%self.id)
        self.stats = obj_from_url(Ipy.API_URL+'/metagenome_statistics/'+self.id+'?verbosity=full', self._auth)

class MetagenomeDisplay(object):
    def __init__(self, mg, def_name=None):
        self.mg = mg
        # hack to get variable name
        if def_name == None:
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self._defined_name = def_name
    
    def annotation(self, annotation='organism', level='domain', source='Subsystems', parent=None):
        if self.mg.stats is None:
            self.mg._set_statistics()
        sub_ann = ''
        if annotation == 'organism':
            annotation = 'taxonomy'
            sub_ann = level
        elif annotation == 'function':
            annotation = 'ontology'
            sub_ann = source
        names  = get_taxonomy(level, parent) if (annotation == 'taxonomy') and (parent is not None) else None
        colors = google_palette(len(self.mg.stats[annotation][sub_ann]))
        data   = []
        for i, d in enumerate(self.mg.stats[annotation][sub_ann]):
            if (names is not None) and (d[0] not in names):
                continue
            data.append({'name': d[0], 'data': [int(d[1])], 'fill': colors[i]})
        annMax  = len(max(self.mg.stats[annotation][sub_ann], key=len))
        pwidth  = 300;
    	pheight = 300;
    	lwidth  = max(300, int(annMax * 7.5));
    	lheight = len(data) * 23;
    	width   = pwidth+lwidth;
    	height  = min(lheight, pheight+(pheight/2)) if lheight > pheight else pheight;
        keyArgs = { 'btype': 'pie',
                    'x_labels': [""],
                    'title': sub_ann,
                    'target': self.mg.id+"_"+level+'_'+random_str(),
                    'show_legend': True,
                    'legendArea': [pwidth+40, 20, lwidth, lheight],
    		        'chartArea': [25, 20, pwidth, pheight],
    		        'width': width,
    		        'height': height,
                    'data': data }
        if annotation == 'taxonomy':
            qname = self._defined_name.replace("'", "\\\'")
            keyArgs['onclick'] = '%s.annotation(annotation="organism", level="%s", parent="\'+params[\'series\']+\'")'%(qname, child_level(level, htype='taxonomy'))
        if Ipy.DEBUG:
            print keyArgs
        else:
            try:
                Ipy.RETINA.graph(**keyArgs)
            except:
                sys.stderr.write("Error producing %s chart"%annotation)
    
    def summary_piechart(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='summary_piechart')
        except:
            sys.stderr.write("Error producing summary piechart\n")
    
    def summary_stats(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='summary_stats')
        except:
            sys.stderr.write("Error producing summary stats\n")
            
    def annotation_piechart(self, annotation='organism', level='domain'):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='annotation_piechart', annotation=annotation, level=level)
        except:
            sys.stderr.write("Error producing annotation piechart\n")
            
    def bp_histogram(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='bp_histogram')
        except:
            sys.stderr.write("Error producing bp histogram\n")
            
    def drisee(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='drisee')
        except:
            sys.stderr.write("Error producing drisee plot\n")
            
    def kmer(self, kmer='abundance'):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='kmer', kmer=kmer)
        except:
            sys.stderr.write("Error producing kmer plot\n")
            
    def rarefaction(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='rarefaction')
        except:
            sys.stderr.write("Error producing rarefaction plot\n")
            
    def rank_abundance(self, level='domain'):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='rank_abundance', level=level)
        except:
            sys.stderr.write("Error producing rank abundance plot\n")
            
    def mixs(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='mixs')
        except:
            sys.stderr.write("Error producing mixs metadata table\n")
            
    def metadata(self):
        try:
            Ipy.RETINA.metagenome(metagenome=self.mg, view='metadata')
        except:
            sys.stderr.write("Error producing full metadata table\n")
