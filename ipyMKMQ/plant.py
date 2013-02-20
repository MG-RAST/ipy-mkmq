#!/usr/bin/env python

import hashlib, traceback
import expression, genopheno, networks, ontology
from ipyTools import *

def get_plant_set(gids=[], def_name=None):
    """Wrapper for Plant object creation, checks if cache (created through unique option set) exists first and returns that.
    returns dict: key = plant_id, value = Plant() object
    
    see: help(Plant)
    """
    if not gids:
        sys.stderr.write("No ids inputted\n")
        return
    cache_id  = "_".join(sorted(gids))
    cache_md5 = hashlib.md5(cache_id).hexdigest()
    cache_obj = load_object(cache_md5)
    if cache_obj is not None:
        print "Loading Plants for selected genomes from cached object"
        return cache_obj
    else:
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        print "Loading Plants for selected genomes through API. Please wait, this may take several minutes ..."
        new_obj = dict([(x, Plant(genome_id=x, def_name="%s['%s']"%(def_name, x))) for x in gids])
        save_object(new_obj, cache_md5)
        print "Done loading through API"
        return new_obj

# class for plants
class Plant(object):
    """Class for working with Plant genomes
    
    self.genome_id   : kbase genome id of plant
    self.EXPRESSION  : expression.PlantExpression object
    self.GENOPHENO   : genopheno.Genotype_PhenotypeAPI object
    self.NETWORKS    : networks.KBaseNetworks object
    self.ONTOLOGY    : ontology.Ontology object
    self.experiments : result of self.GENOPHENO.get_experiments(self.genome_id)
    self.traits      : result of self.GENOPHENO.get_traits(self.experiments[0])
    
    """
    def __init__(self, genome_id, def_name=None):
        self.genome_id   = genome_id
        self.EXPRESSION  = expression.PlantExpression(Ipy.EXPRESSION_URL)
        self.GENOPHENO   = genopheno.Genotype_PhenotypeAPI(Ipy.GENOPHENO_URL)
        self.NETWORKS    = networks.KBaseNetworks(Ipy.NETWORKS_URL)
        self.ONTOLOGY    = ontology.Ontology(Ipy.ONTOLOGY_URL)
        self.experiments = self.GENOPHENO.get_experiments(self.genome_id)
        self.traits      = self.GENOPHENO.get_traits(self.experiments[0])
        # hack to get variable name
        if def_name == None:
            (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
            def_name = text[:text.find('=')].strip()
        self.defined_name = def_name

    def get_variations(self, count=5):
        return self.GENOPHENO.traits_to_variations(self.traits[0], count)

    def plot_variations(self, count=5, variations=None, title='', width=1100, height=400, x_min=0, x_max=None, y_min=0, y_max=10, arg_list=False):
        if not variations:
            variations = self.get_variations(count)
        colors  = google_palette(count)
        series  = []
        points  = []
        lengths = []
        offsets = []
        if not title:
            title = "Manhattan Plot for %s"%variations["trait"]["trait_name"]
        for i in range(count):
            series.append({ "name": str(i+1), "color": colors[i], "shape": "circle"})
            points.append([])
            lengths.append(0)
            offsets.append(0)
        for i in variations["variations"]:
            lengths[i[0]] = max(lengths[i[0]], i[1])
        for i in range(len(lengths)):
            if (i == 0):
                offsets[i] = 10000
            else:
                offsets[i] = offsets[i-1] + lengths[i-1] + 1000000
        for i in variations["variations"]:
            points[i[0]].append({ "x": i[1] + offsets[i[0]], "y": i[2] })
        if not x_max:
            x_max = offsets[count-1] + lengths[count-1]
        keyArgs = { 'width': width,
                    'height': height,
                    'x_min': x_min,
                    'x_max': x_max,
                    'y_min': y_min,
                    'y_max': y_max,
                    'connected': False,
                    'show_dots': True,
                    'data': {"series": series, "points": points}
                  }
        if Ipy.DEBUG:
            print keyArgs
        if arg_list:
            return keyArgs
        else:
            try:
                Ipy.RETINA.plot(**keyArgs)
            except:
                sys.stderr.write("Error producing manhattan plot\n")
            return None
