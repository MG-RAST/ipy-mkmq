#!/usr/bin/env python

from ipyTools import *

class Metagenome:
    def __init__(self, mgid, metadata=True, stats=False):
        # get mg
        metagenome = self._get_metagenome(metadata)
        if metagenome is not None:
            for key, val in metagenome.iteritems():
                self.__setitem__(key, val)
        else:
            self.id = mgid
            self.name = None
        # get stats
        self.stats = None
        if stats:
            self.stats = self._get_statistics()
        
    def _get_metagenome(self, metadata):
        verb = 'full' if metadata else 'verbose'
        return obj_from_url(API_URL+'metagenome/'+self.mgid+'?verbosity='+verb)

    def _get_statistics(self):
        return obj_from_url(API_URL+'metagenome_statistics/'+self.mgid)