#!/usr/bin/env python

import os, shutil, traceback
from collection import Collection
from ipyTools import *

class Project(Collection):
    """Class representation of Project object:
         "id"             : [ 'string', 'unique object identifier' ],
         "name"           : [ 'string', 'human readable identifier' ],
         "libraries"      : [ 'list',  ['reference library', 'a list of references to the related library objects'] ],
         "samples"        : [ 'list',  ['reference sample', 'a list of references to the related sample objects'] ],
         "analyzed"       : [ 'list',  ['reference metagenome', 'a list of references to the related metagenome objects'] ],
         "description"    : [ 'string', 'a short, comprehensive description of the project' ],
         "funding_source" : [ 'string', 'the official name of the source of funding of this project' ],
         "pi"             : [ 'string', 'the first and last name of the principal investigator of the project' ],
         "metadata"       : [ 'hash',   'key value pairs describing metadata' ],
         "created"        : [ 'date',   'time the object was first created' ],
         "version"        : [ 'integer','version of the object' ],
         "url"            : [ 'uri',    'resource location of this object instance' ],
         "status"         : [ 'cv',     [ ['public', 'object is public'],
        						           ['private', 'object is private'] ] ]
    """
    def __init__(self, pid, metadata=True, stats=True, auth=None, def_name=None, cache=False, reset_cache=False):
        # set project
        self.cache = Ipy.NB_DIR+'/'+pid if cache else None
        project = None
        # hack to get variable name
        if def_name == None:
            try:
                (filename,line_number,function_name,text)=traceback.extract_stack()[-2]
                def_name = text[:text.find('=')].strip()
            except:
                pass
        self.defined_name = def_name        
        # reset cache if asked
        if reset_cache and os.path.isdir(self.cache):
            shutil.rmtree(self.cache)
        # make cache dir
        if self.cache and (not os.path.isdir(self.cache)):
            os.mkdir(self.cache)
        # try load from cached
        if self.cache and os.path.isdir(self.cache) and os.path.isfile(self.cache+'/'+pid+'.json'):
            try:
                project = json.load(open(self.cache+'/'+pid+'.json', 'rU'))
                sys.stdout.write("project '%s' loaded from cache %s\n"%(self.defined_name, pid))
            except:
                pass
        # load from api
        if project is None:
            project = self._get_project(pid, metadata, auth)
            if project and self.cache and os.path.isdir(self.cache):
                # cache it if dir given and not loaded from file
                try:
                    json.dump(project, open(self.cache+'/'+pid+'.json', 'w'))
                    sys.stdout.write("project '%s' saved to cache %s\n"%(self.defined_name, pid))
                except:
                    pass
        if project is None:
            self.id = pid
            self.name = None
            return
        for key, val in project.iteritems():
            setattr(self, key, val)
        # call collection init - from cache if given
        Collection.__init__(self, self.mgids(), metadata=metadata, stats=stats, auth=auth, def_name=self.defined_name, cache=self.cache)
    
    def _get_project(self, pid, metadata, auth):
        verb = 'full' if metadata else 'verbose'
        return obj_from_url(Ipy.API_URL+'project/'+pid+'?verbosity='+verb, auth)

    def mgids(self):
        mlist = []
        if hasattr(self, 'analyzed'):
            mlist = map(lambda x: x[0], self.analyzed)
        return mlist
