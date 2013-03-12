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
        _mgids : [ 'list', 'inputted metagenome ids' ]
    
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
        self._defined_name = def_name
        # get metagenomes
        self.metagenomes = self._get_metagenomes(cache)
        self.rarefaction = Rarefaction(mgObjs=self.metagenomes.values())
    
    def _get_metagenomes(self, cache):
        mgs = {}
        for mg in self._mgids:
            keyArgs = { 'stats': self._stats,
                        'auth': self._auth,
                        'cache': cache,
                        'def_name': '%s.metagenomes["%s"]'%(self._defined_name, mg)
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
    
    def sub_mgs(self, category=None, field=None, value=None):
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
    
    def show_metadata(self, mgids=None, arg_list=False):
        header = []
        tdata  = []
        mdata  = dict([(x, {}) for x in Ipy.MD_CATS])
        submgs = mgids if mgids and (len(mgids) > 0) else self.mgids()
        for mid in submgs:
            if (mid in self.metagenomes) and hasattr(self.metagenomes[mid], 'metadata'):
                header.append(mid)
        if len(header) == 0:
            sys.stderr.write("No metadata to display\n")
        for i, mid in enumerate(header):
            for cat, data in self.metagenomes[mid].metadata.iteritems():
                for field, value in data['data'].iteritems():
                    if field not in mdata[cat]:
                        mdata[cat][field] = ['' for x in range(len(header))]
                    mdata[cat][field][i] = value
        for cat in mdata.iterkeys():
            for field in mdata[cat].iterkeys():
                tdata.append( [cat, field] + mdata[cat][field] )
        keyArgs = { 'width': 700,
                    'height': 600,
                    'target': "metadata_table_"+random_str(),
                    'data': {'data': tdata, 'header': ['category', 'field'] + header},
                    'rows_per_page': 20 }
        if Ipy.DEBUG:
            print keyArgs
        if arg_list:
            return keyArgs
        else:
            try:
                Ipy.RETINA.table(**keyArgs)
            except:
                sys.stderr.write("Error producing metadata table\n")
            return None

    def plot_rarefaction(self, mgids=None, width=800, height=300, title="", x_title="", y_title="", legend=True, arg_list=False):
        return self.rarefaction.plot(mgids=mgids, width=width, height=height, title=title, x_title=x_title, y_title=y_title, legend=legend, arg_list=arg_list)

    def barchart_taxon(self, level='domain', parent=None, width=800, height=0, x_rotate='0', title="", legend=True):
        children = get_taxonomy(level, parent) if parent is not None else None
        self._barchart('taxonomy', level, width, height, x_rotate, title, legend, names=children)

    def barchart_function(self, source='Subsystems', width=800, height=0, x_rotate='0', title="", legend=True):
        self._barchart('ontology', source, width, height, x_rotate, title, legend)

    def _barchart(self, atype, level, width, height, x_rotate, title, legend, names=None):
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
        height  = height if height else len(annL)*len(self.metagenomes)*7.5
        lheight = min(height, len(self.metagenomes)*35)
        lwidth  = len(max(annL, key=len)) * 7.2
        cwidth  = 0.85 if legend else 0.99
        keyArgs = { 'btype': 'row',
                    'width': width+lwidth,
                    'height': height,
                    'x_labels': annL,
                    'x_labels_rotation': x_rotate,
                    'title': title,
                    'target': '_'.join(self.mgids())+"_"+level+'_'+random_str(),
                    'show_legend': legend,
                    'legendArea': [0.87, 0.05, 0.2, lheight],
                    'chartArea': [lwidth, 0.02, cwidth, 0.95],
                    'data': data }
        if Ipy.DEBUG:
            print keyArgs
        try:
            Ipy.RETINA.graph(**keyArgs)
        except:
            sys.stderr.write("Error producing chart")
            return None
