import os.path as op

from transform import EBSTransform, TransformPoints
from shape import FindSpots,Thinning, SBR
from template import TemplateSource
from phantom import PhantomGenerate

_root = op.abspath( op.join( op.dirname(__file__), '..', '..' ) )

def get_datapath():
    return op.join( _root, 'data' )
