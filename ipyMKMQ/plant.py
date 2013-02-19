#!/usr/bin/env python

import expression, genopheno, networks, ontology
from ipyTools import *
#import cdmi

# class for plants
class Plant(object):
    def __init__(self, genome_id):
        self.EXPRESSION  = expression.PlantExpression(Ipy.EXPRESSION_URL)
        self.GENOPHENO   = genopheno.Genotype_PhenotypeAPI(Ipy.GENOPHENO_URL)
        self.NETWORKS    = networks.KBaseNetworks(Ipy.NETWORKS_URL)
        self.ONTOLOGY    = ontology.Ontology(Ipy.ONTOLOGY_URL)
        self.experiments = self.GENOPHENO.get_experiments(genome_id)
        self.traits      = self.GENOPHENO.get_traits(self.experiments[0])

    def get_variations(self, num):
        return self.GENOPHENO.traits_to_variations(self.traits[0], num)

    def plot_variations(self, num, title='', width=1100, height=400, x_min=0, x_max=None, y_min=0, y_max=10, arg_list=False):
        variations = self.get_variations(num)
        colors  = google_palette(num)
        series  = []
        points  = []
        lengths = []
        offsets = []
        if not title:
            title = "Manhattan Plot for %s"%variations["trait"]["trait_name"]
        for i in range(num):
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
            x_max = offsets[num-1] + lengths[num-1]
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
