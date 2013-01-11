#!/usr/bin/env python

from collection import Collection
from ipyTools import *

class Project(Collection):
    def __init__(self, pid, metadata=True, stats=True, auth=None, def_name=None):
        # set project
        project = self._get_project(pid, metadata, auth)
        if project is None:
            self.id = pid
            self.name = None
            return
        for key, val in project.iteritems():
            setattr(self, key, val)
        # call collection init
        super(Collection, self).__init__(self.metagenomes(), metadata, stats, auth, def_name)
    
    def _get_project(self, pid, metadata, auth):
        verb = 'full' if metadata else 'verbose'
        key  = '&auth='+auth if auth else ''
        return obj_from_url(Ipy.API_URL+'project/'+pid+'?verbosity='+verb+key)
    
    def metagenomes(self):
        mlist = []
        if hasattr(self, 'analyzed'):
            mlist = map(lambda x: x[0], self.analyzed)
        return mlist
