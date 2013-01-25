"""
prototype ipy tools.

"""
__author__ = 'Travis Harrison'
__version__ = '0.2'
__description__ = 'iPython Tools for Qiime-Matr-QC'
__all__ = ["analysis","qc","ipyTools","flotplot","retina","metagenome","project","collection"]

from analysis import Analysis, AnalysisSet
from project import Project
from collection import Collection
from metagenome import Metagenome
from qc import QC, Drisee, NucleoProfile, Rarefaction, merge_drisee_profile
from retina import Retina
from flotplot import FlotPlot