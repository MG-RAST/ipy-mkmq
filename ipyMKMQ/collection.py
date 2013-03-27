#!/usr/bin/env python

import sys, os, hashlib, traceback
from collections import defaultdict
from metagenome import Metagenome
from ipyTools import *
from qc import Rarefaction

def get_collection(mgids=[], auth=None, stats=True, def_name=None):
    """Wrapper for Collection object creation, checks if cache (created through unique option set) exists first and returns that.
    
    see: help(Collection)
    """
    if not mgids:
        sys.stderr.write("No ids inputted\n")
        return
    cache_id  = "_".join(sorted(mgids))+"_"+('1' if stats else '0')
    cache_md5 = hashlib.md5(cache_id).hexdigest()
    cache_obj = load_object(cache_md5)
    if cache_obj is not None:
        print "Loading Collection for selected metagenomes from cached object"
        return cache_obj
    else:
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        print "Loading Collection for selected metagenomes through API. Please wait, this may take several minutes ..."
        new_obj = Collection(mgids=mgids, auth=auth, stats=stats, def_name=def_name)
        save_object(new_obj, cache_md5)
        print "Done loading through API"
        return new_obj

class Collection(object):
    """Class representation of Collection object:
        metagenomes : [ 'hash', 'key = metagenome_id, value = Metagenome() object']
        rarefaction : Rarefaction object for collection metagenomes
        _mgids      : [ 'list', 'inputted metagenome ids' ]
        display     : 'CollectionDisplay Object - help(this_name.display)'
        
    see: help(Metagenome)
    """
    def __init__(self, mgids=[], stats=True, auth=None, def_name=None, cache=False):
        self._auth  = auth
        self._stats = stats
        self._mgids = mgids
        # hack to get variable name
        if def_name == None:
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self.defined_name = def_name
        # get metagenomes
        self.metagenomes = self._get_metagenomes(cache)
        # set display
        self.display = CollectionDisplay(self.metagenomes.values(), self.defined_name+'.display')
    
    def _get_metagenomes(self, cache):
        mgs = {}
        for mg in self._mgids:
            keyArgs = { 'stats': self._stats,
                        'auth': self._auth,
                        'cache': cache,
                        'def_name': '%s.metagenomes["%s"]'%(self.defined_name, mg)
                       }
            mgs[mg] = Metagenome(mg, **keyArgs)
        return mgs
    
    def _set_statistics(self):
        self._stats = True
        for mg in self.metagenomes.itervalues():
            mg._set_statistics()
    
    def mgids(self):
        return self._mgids
    
    def get_stat(self, mgid=None, stat=None, mgid_set=[]):
        if not (mgid and stat and (mgid in self._mgids)):
            return []
        if not self._stats:
            self._set_statistics()
        if stat not in self.metagenomes[mgid].stats['sequence_stats']:
            return []
        stat_list = [ toNum(self.metagenomes[mgid].stats['sequence_stats'][stat]) ]
        if not mgid_set:
            mgid_set = self._mgids
        for m in mgid_set:
            if m == mgid:
                continue
            if stat in self.metagenomes[m].stats['sequence_stats']:
                stat_list.append( toNum(self.metagenomes[m].stats['sequence_stats'][stat]) )
        return stat_list

    def metadata_fields(self, table=True):
            tdata = []
            mdata = dict([(x, set()) for x in Ipy.MD_CATS])
            for mg in self.metagenomes.itervalues():
                if not hasattr(mg, 'metadata'):
                    continue
                for cat in Ipy.MD_CATS:
                    if cat not in mg.metadata:
                        continue
                    for key in mg.metadata[cat]['data'].iterkeys():
                        mdata[cat].add(key)
            if not table:
                return mdata
            for cat in mdata.iterkeys():
                for field in sorted(mdata[cat]):
                    tdata.append([cat, field])
            keyArgs = { 'width': 400,
                        'height': 600,
                        'target': 'metadata_fields_'+random_str(),
                        'data': {'data': tdata, 'header': ['category', 'field']},
                        'rows_per_page': 20 }
            if Ipy.DEBUG:
                print keyArgs
            try:
                Ipy.RETINA.table(**keyArgs)
            except:
                sys.stderr.write("Error producing metadata table\n")
    
    def search_metadata(self, category=None, field=None, value=None):
        sub_mgs = set()
        all_fields = []
        for f in self.metadata_fields(table=False).itervalues():
            all_fields.extend(list(f))
        if not (category and (category in Ipy.MD_CATS)):
            sys.stderr.write("category must be one of: %s\n"%", ".join(Ipy.MD_CATS))
            return self.mgids()
        if not (field and value and (field in all_fields)):
            sys.stderr.write("field '%s' does not exist\n"%field)
            return self.mgids()
        for mid, mg in self.metagenomes.iteritems():
            if not (hasattr(mg, 'metadata') and (category in mg.metadata)):
                continue
            for key, val in mg.metadata[category]['data'].iteritems():
                if key == field:
                    x = str(val).find(value)
                    if x != -1:
                        sub_mgs.add(mid)
        return list(sub_mgs)

class CollectionDisplay(object):
    """Class containing functions to display metagenome collection visualizations:
        annotation        : interactive barchart of organism or functional abundances with clickable drilldown
        annotation_chart  : static barchart of organism or functional abundances
        drisee            : plot of DRISEE sequencing error
        kmer              : plot of kmer profile
        metadata          : interactive table of full metadata
        mixs              : table of GSC MIxS metadata
        rarefaction       : plot of species rarefaction
        summary_chart     : barchart of summary sequence hits
        summary_stats     : table of summary statistics
    """
    def __init__(self, mgs, def_name=None):
        self.mgs = mgs
        self.mgids = map(lambda x: x.id, mgs)
        # hack to get variable name
        if def_name == None:
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self.defined_name = def_name

    def annotation(self, annotation='organism', level='domain', source='Subsystems', title='', parent=None, arg_list=False):
        sub_ann = ''
        if annotation == 'organism':
            annotation = 'taxonomy'
            sub_ann = level
        elif annotation == 'function':
            annotation = 'ontology'
            sub_ann = source
        names  = get_taxonomy(level, parent) if (annotation == 'taxonomy') and (parent is not None) else None
        colors = google_palette(len(self.mgs))
        data = []
        annD = {}
        for i, mg in enumerate(self.mgs):
            data.append({'name': mg.id, 'data': [], 'fill': colors[i]})
            for d in mg.stats[annotation][level]:
                if (names is not None) and (d[0] not in names):
                    continue
                annD[ d[0] ] = 1
        annL = sorted(annD.keys())
        for i, d in enumerate(data):
            annMG = {}
            for a, v in self.mgs[i].stats[annotation][level]:
                annMG[a] = v
            for a in annL:
                if a in annMG:
                    d['data'].append(int(annMG[a]))
                else:
                    d['data'].append(0)
        
        width   = 600
        height  = len(annL) * len(self.mgs) * 7.5
        lwidth  = len(max(annL, key=len)) * 7.8
        lheight = len(self.mgs) * 35
        if height < 100:
            height = 100
        keyArgs = { 'title': sub_ann if not title else title,
        	        'btype': 'row',
        		    'x_labels': annL,
        		    'target': '_'.join(self.mgids)+"_"+level+'_'+random_str(),
        		    'show_legend': True,
        		    'legendArea': [width+lwidth, 50, 150, lheight],
        		    'chartArea': [lwidth, 50, 0.81, height],
        		    'width': width+lwidth+150,
        		    'height': max(height, lheight),
        		    'data': data }
        if Ipy.DEBUG:
            print keyArgs
        if arg_list:
            return keyArgs
        else:
            try:
                Ipy.RETINA.graph(**keyArgs)
            except:
                sys.stderr.write("Error producing %s chart"%annotation)
            return None
        
    def summary_chart(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='summary_chart', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing summary chart\n")

    def summary_stats(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='summary_stats', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing summary stats\n")
            
    def annotation_chart(self, annotation='organism', level='domain', arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='annotation_chart', annotation=annotation, level=level, arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing annotation chart\n")
            
    def drisee(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='drisee', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing drisee plot\n")
            
    def kmer(self, kmer='abundance', arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='kmer', kmer=kmer, arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing kmer plot\n")
            
    def rarefaction(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='rarefaction', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing rarefaction plot\n")
            
    def mixs(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='mixs', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing mixs metadata table\n")
            
    def metadata(self, arg_list=False, target=None):
        try:
            Ipy.RETINA.collection(metagenomes=self.mgs, view='metadata', arg_list=arg_list, target=target)
        except:
            sys.stderr.write("Error producing full metadata table\n")
