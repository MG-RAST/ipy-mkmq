#!/usr/bin/env python

import expression, genopheno, networks, ontology
import config
#import cdmi

# class for plants
class Plant(object):
    """Constants for plant interface"""
    EXPRESSION = expression.PlantExpression(config.EXPRESSION_URL)
    GENOPHENO  = genopheno.Genotype_PhenotypeAPI(config.GENOPHENO_URL)
    NETWORKS   = networks.KBaseNetworks(config.NETWORKS_URL)
    ONTOLOGY   = ontology.Ontology(config.ONTOLOGY_URL)
