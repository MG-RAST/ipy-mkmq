#!/usr/bin/env python

from ipyTools import *

class Metagenome:
    def __init__(self, mgid):
        metagenome = self.get_metagenome()
        if metagenome not None:
            for key, val in metagenome.iteritems():
                self.__setitem__(key, val)
        else:
            self.id = mgid
            self.name = None
        
    def get_metagenome(self):
        return obj_from_url(API_URL+'metagenome/'+self.mgid)