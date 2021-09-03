import imp
import os
import inspect
import bpy
import bmesh

from pro import context

from .material import MaterialManager

from .util import VertexRegistry

from .op_decompose import Decompose
from .op_split import Split
from .op_extrude import Extrude
from .op_extrude2 import Extrude2
from .op_color import Color
from .op_material import Material
from .op_texture import Texture
from .op_delete import Delete
from .op_join import Join
from .op_inset import Inset
from .op_inset2 import Inset2
from .op_rectangle import Rectangle
from .op_hip_roof import HipRoof
from .op_complex_hip_roof import ComplexHipRoof
from .op_dutch_roof import DutchRoof
from .op_gable_roof import GableRoof
from .op_hip_roof2 import HipRoof2
from .op_copy import Copy
from .op_translate import Translate



from pro.base import Param

from .shape import getInitialShape

from .join import JoinManager


def buildFactory():
    factory = context.factory
    factory["Decompose"] = Decompose
    factory["Split"] = Split
    factory["Extrude"] = Extrude
    factory["Extrude2"] = Extrude2
    factory["Color"] = Color
    factory["Material"] = Material
    factory["Texture"] = Texture
    factory["Delete"] = Delete
    factory["Join"] = Join
    factory["Inset"] = Inset
    factory["Inset2"] = Inset2
    factory["Rectangle"] = Rectangle
    factory["HipRoof"] = HipRoof
    factory["ComplexHipRoof"]= ComplexHipRoof
    factory["DutchRoof"]= DutchRoof
    factory["GableRoof"]= GableRoof
    factory["HipRoof2"]= HipRoof2
    factory["Copy"] = Copy
    factory["Translate"] = Translate


def apply(ruleFile, startRule="Begin"):
    from .bl_util import create_rectangle

    blenderContext = context.blenderContext
    obj = blenderContext.object
    if obj:
        bpy.ops.object.mode_set(mode="OBJECT")

    noMeshCondition = not obj or obj.type != "MESH"
    if noMeshCondition or len(obj.data.polygons) != 1:
        if not noMeshCondition:
            # delete if it's non-flat mesh
            bpy.ops.object.delete()
        create_rectangle(blenderContext, 20, 10)
    # apply all transformations to the active Blender object
    bpy.ops.object.transform_apply(location=True, rotation=True, scale=True)
    # setting the path to the rule for context
    context.ruleFile = ruleFile if isinstance(
        ruleFile, str) else ruleFile.__file__
    params = None
    mesh = blenderContext.object.data
    # create a uv layer for the mesh
    # TODO: a separate pass through the rules is needed to find out how many uv layers are needed
    mesh.uv_layers.new(name=Texture.defaultLayer)
    # initialize the context
    context.init()
    # initializing bmesh instance
    bm = bmesh.new()
    bm.from_mesh(mesh)
    if hasattr(bm.faces, "ensure_lookup_table"):
        bm.faces.ensure_lookup_table()
    context.addAttribute("bm", bm)
    # list of unused faces for removal
    context.addAttribute("facesForRemoval", [])
    # set up the material registry
    context.addAttribute("materialManager", MaterialManager())
    # set up vertex registry to ensure vertex uniqueness
    context.addAttribute("vertexRegistry", VertexRegistry())
    # set a constructor for join manager, it may be replaced by actual instance of the join manager
    context.addAttribute("joinManager", JoinManager)

    # push the initial state with the initial shape to the execution stack
    context.pushState(shape=getInitialShape(bm))

    if isinstance(ruleFile, str):
        module = getModule(ruleFile)

        # prepare context internal stuff
        context.prepare()
        # params is a list of tuples: (paramName, instanceofParamClass)
        params = getParams(module)
    else:
        # ruleFile is actually a module
        module = ruleFile

    # setting the current operator to a dummy one to avoid an exception
    class dummy:
        def addChildOperator(self, o): pass

        def removeChildOperators(self, numParts): pass
    context.operator = dummy()
    # evaluate the rule set: this command actually read and parse the rule set --> create 3D buildings according to rules 
  
    #execute first rule (that calls the other rules). For instance for test.py startrule Begin is called and run the whole rulefile 
    getattr(module, startRule)().execute()
    # remove unused faces from context.facesForRemoval
    bmesh.ops.delete(bm, geom=context.facesForRemoval, context='FACES')
    # there still may be some doubles, inspite of the use of util.VertexMaterial
    #bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.0001)

    # clean up context.facesForRemoval
    context.facesForRemoval = []
    context.executeDeferred()
    # remove unused faces from context.facesForRemoval
    bmesh.ops.delete(bm, geom=context.facesForRemoval, context='FACES')

    # write everything back to the mesh
    bm.to_mesh(mesh)
    
    # for each param we create a custom property to store the parameter value for the active object (building)
    # Now you can retrieve any parameter, i.e. bpy.data.objects['BCGA']['height'] gives the building height
    for param in params:
        param_name=param[0]
        param_value=param[1]
        obj[param_name]= param_value.getValue()
       
    # cleaning context from blender specific members
    context.removeAttributes()

    return (module, params)


def isParam(member):
    """A predicate for the inspect.getmembers call"""
    return isinstance(member, Param)


def getModule(ruleFile):
    """Returns Python module object given a path to the rule file"""
    # remove extension from ruleFile if it was provided
    ruleFile = os.path.splitext(ruleFile)[0]
    moduleName = os.path.basename(ruleFile)
    _file, _pathname, _description = imp.find_module(
        moduleName, [os.path.dirname(ruleFile)])
    module = imp.load_module(moduleName, _file, _pathname, _description)
    _file.close()
    return module


def getParams(module):
    """Returns a list of tuples: (paramName, instanceofParamClass)"""
    return [m for m in inspect.getmembers(module, isParam)]


buildFactory()
