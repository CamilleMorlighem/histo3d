import math
import os
import sys

bl_info = {
	"name": "BCGA",
	"author": "Vladimir Elistratov <vladimir.elistratov@gmail.com>",
	"version": (1, 0, 0),
	"blender": (2, 80, 0),
	"location": "View3D > Tool Shelf",
	"description": "BCGA: Computer Generated Architecture for Blender",
	"warning": "",
	"wiki_url": "https://github.com/vvoovv/bcga/wiki",
	"tracker_url": "https://github.com/vvoovv/bcga/issues",
	"support": "COMMUNITY",
	"category": "BCGA",
}

numFloatParams = 200
numColorParams = 50

for path in sys.path:
	if "bpro" in path:
		path = None
		break
if path:
	# we need to add path to bpro package
	sys.path.append(os.path.dirname(__file__))

import bpy
import bpro

from pro import context as proContext
from pro.base import ParamFloat, ParamColor

from bpro.bl_util import create_rectangle, align_view, first_edge_ymin

#### Returns full path to a BDGA script or None if it does not exist #### 
def getRuleFile(ruleFile, operator):
	"""
	Returns full path to a BCGA script or None if it does not exist.
	"""
	# Checks if the path specified looks like a path (dans la forme)
	if len(ruleFile)>1 and ruleFile[:2]=="//":
		ruleFile = ruleFile[2:]
	ruleFile = os.path.join(os.path.dirname(bpy.data.filepath), ruleFile)
	# Checks if the path exists in the computer 
	if not os.path.isfile(ruleFile):
		operator.report({"ERROR"}, "The BCGA script '%s' not found" % ruleFile)
		ruleFile = None
	return ruleFile

# Add custom properties (string type) to the scene (all scenes now have those properties). note that they should appear in 
# scene custom properties in the UI but it does not show in blender 2.83 in the custom properties pannel (known issue) 
bpy.types.Scene.bcgaScript = bpy.props.StringProperty(
	name = "Script",
	description = "Path to a BCGA script",
	subtype="FILE_PATH"
)

bpy.types.Scene.bakingBcgaScript = bpy.props.StringProperty(
	name = "Low poly script",
	description = "Path to a BCGA script with a low poly model",
	subtype="FILE_PATH"
)

#?
class CustomFloatProperty(bpy.types.PropertyGroup):
	"""A bpy.types.PropertyGroup descendant for bpy.props.CollectionProperty"""
	value: bpy.props.FloatProperty(name="")
#?
class CustomColorProperty(bpy.types.PropertyGroup):
	"""A bpy.types.PropertyGroup descendant for bpy.props.CollectionProperty"""
	value: bpy.props.FloatVectorProperty(name="", subtype='COLOR', min=0.0, max=1.0)


#### Draw the main BCGA panel in the UI #### 
class ProMainPanel(bpy.types.Panel):
	bl_label = "Main"
	bl_space_type = "VIEW_3D"
	bl_region_type = "UI"
	#bl_context = "objectmode"
	bl_category = "BCGA"
	
	def draw(self, context):
		scene = context.scene
		layout = self.layout
		# Add a button in the panel to run the "object.footprint_set" operator 
		layout.row().operator_menu_enum("object.footprint_set", "size", text="Footprint")
		layout.separator()
		# Add a button in the panel to set the custom propery "bcgaScript" (we actually searh for that property)
		layout.row().prop_search(scene, "bcgaScript", bpy.data, "texts")
		# Add a button in the panel to run the "object.apply_pro_script" operator 
		layout.row().operator("object.apply_pro_script")


#### Draw the baking panel in the UI #### 
class BakingPanel(bpy.types.Panel):
	bl_label = "Baking"
	bl_space_type = "VIEW_3D"
	bl_region_type = "UI"
	bl_category = "BCGA"
	#bl_options = {"DEFAULT_CLOSED"}
	
	def draw(self, context):
		scene = context.scene
		layout = self.layout
		# Add a button in the panel to set the custom propery "bakingBcgaScript" (we actually searh for that property)
		layout.row().prop_search(scene, "bakingBcgaScript", bpy.data, "texts")
		# Add a button in the panel to run the "object.bake_pro_model" operator 
		self.layout.operator("object.bake_pro_model")


#### Draw the first edge panel in the UI #### 
class FirstEdgePanel(bpy.types.Panel):
	bl_label = "First edge"
	bl_space_type = "VIEW_3D"
	bl_region_type = "UI"
	bl_category = "BCGA"
	bl_options = {"DEFAULT_CLOSED"}
	
	def draw(self, context):
		# Add a button in the panel to run the "object.first_edge_ymin" operator
		self.layout.operator("object.first_edge_ymin")


#### Create operator that creates 3D building from rule file and footprint #### 
class Pro(bpy.types.Operator):
	bl_idname = "object.apply_pro_script"
	bl_label = "Apply"
	bl_options = {"REGISTER", "UNDO"}
	
	collectionFloat: bpy.props.CollectionProperty(type=CustomFloatProperty)
	collectionColor: bpy.props.CollectionProperty(type=CustomColorProperty)
	
	initialized = False
	
	def initialize(self):
		if self.initialized:
			return
		for _ in range(numFloatParams):
			self.collectionFloat.add()
		for _ in range(numColorParams):
			self.collectionColor.add()
		self.initialized = True
	
	def invoke(self, context, event):
		self.initialize()
		proContext.blenderContext = context
		#ruleFile = bpy.data.texts[context.scene.bcgaScript].filepath
		ruleFile = context.scene.bcgaScript
		if ruleFile:
			# append the directory of the ruleFile to sys.path
			ruleFileDirectory = os.path.dirname(os.path.realpath(os.path.expanduser(ruleFile)))
			if ruleFileDirectory not in sys.path:
				sys.path.append(ruleFileDirectory)

			# retrieve the parameters from the rulefile (height,...)
			# calls bpro.apply which runs the ruleFile 
			module,params = bpro.apply(ruleFile)
			
			#align_view(context.object)
			
			self.module = module
			self.params = params

	
			numFloats = 0
			numColors = 0
			# for each entry in self.params create a new item in self.collection
			# the following displays the parameters into the redo pannel when u call the operator 
			for param in self.params:
				param = param[1]
				if isinstance(param, ParamFloat):
					collectionItem = self.collectionFloat[numFloats]
					numFloats += 1
				elif isinstance(param, ParamColor):
					collectionItem = self.collectionColor[numColors]
					numColors += 1
				collectionItem.value = param.getValue()
				param.collectionItem = collectionItem
	
		return {"FINISHED"}
	
	def execute(self, context):
		self.initialize()
		proContext.blenderContext = context
		#ruleFile = bpy.data.texts[context.scene.bcgaScript].filepath
		ruleFile = context.scene.bcgaScript
		if ruleFile:
			# append the directory of the ruleFile to sys.path
			ruleFileDirectory = os.path.dirname(os.path.realpath(os.path.expanduser(ruleFile)))
			if ruleFileDirectory not in sys.path:
				sys.path.append(ruleFileDirectory)
			# retrieve the parameters from the rulefile (height,...)
			try: 
				module,params = bpro.apply(ruleFile)
			
				#align_view(context.object)
				
				self.module = module
				self.params = params
			
				numFloats = 0
				numColors = 0
				# for each entry in self.params create a new item in self.collection
				for param in self.params:
					param = param[1]
					if isinstance(param, ParamFloat):
						collectionItem = self.collectionFloat[numFloats]
						numFloats += 1
					elif isinstance(param, ParamColor):
						collectionItem = self.collectionColor[numColors]
						numColors += 1
					collectionItem.value = param.getValue()
					param.collectionItem = collectionItem
			except: 
				pass
			
			
		return {"FINISHED"}

		"""
		#the following was the content of the execute function before but it raised an error so i removed it and replace 
		#by content of invoke function 
		proContext.blenderContext = context
		for param in self.params:
			param = param[1]
			param.setValue(getattr(param.collectionItem, "value"))
		bpro.apply(self.module)
		
		#align_view(context.object)
		
		return {"FINISHED"}"""
	
	# draw the panel with the parameters and their values so that you can modify the values in the UI 
	def draw(self, context):
		layout = self.layout
		if hasattr(self, "params"):
			# self.params is a list of tuples: (paramName, instanceofParamClass)
			for param in self.params:
				paramName = param[0]
				row = layout.split()
				row.label(text=paramName + ":")
				row.prop(param[1].collectionItem, "value")


class Bake(bpy.types.Operator):
	bl_idname = "object.bake_pro_model"
	bl_label = "Bake"
	bl_options = {"REGISTER", "UNDO"}

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == "CYCLES"

	def execute(self, context):
		proContext.blenderContext = context
		bpy.ops.object.select_all(action="DESELECT")
		# remember the original object, it will be used for low poly model
		lowPolyObject = context.object
		lowPolyObject.select_set(True)
	
		bpy.ops.object.duplicate()
		highPolyObject = context.object

		# high poly model
		#######ruleFile = bpy.data.texts[context.scene.bcgaScript].filepath
		ruleFile = context.scene.bcgaScript

		if ruleFile:
			highPolyParams = bpro.apply(ruleFile)[1]
			# convert highPolyParams to a dict paramName->instanceofParamClass
			highPolyParams = dict(highPolyParams)
			
			# low poly model
			context.view_layer.objects.active = lowPolyObject
			######ruleFile = bpy.data.texts[context.scene.bakingBcgaScript].filepath
			ruleFile = context.scene.bakingBcgaScript
			if ruleFile:
				name = lowPolyObject.name
				module = bpro.getModule(ruleFile)
				lowPolyParams = bpro.getParams(module)
				# Apply highPolyParams to lowPolyParams
				# Normally lowPolyParams is a subset of highPolyParams
				for paramName,param in lowPolyParams:
					if paramName in highPolyParams:
						param.setValue(highPolyParams[paramName].getValue())
				bpro.apply(module)
				# unwrap the low poly model
				bpy.ops.object.mode_set(mode="EDIT")
				bpy.ops.mesh.select_all(action="SELECT")
				bpy.ops.uv.smart_project()
				# prepare settings for baking
				bpy.ops.object.mode_set(mode="OBJECT")
				#highPolyObject.select_set(True)
				bpy.data.objects['BCGA'].select_set(True) #added by cam 
				highPolyObject=bpy.data.objects['BCGA'] #added by cam 
				bpy.context.scene.cycles.bake_type = "DIFFUSE"
				lowPolyObject=bpy.data.objects['BCGA.002'] #added by cam 
				context.view_layer.objects.active = lowPolyObject
				bpy.context.scene.cycles.use_bake_selected_to_active = True
				# create a new image with default settings for baking
				image = bpy.data.images.new(lowPolyObject.name, width=512, height=512, alpha=True)

				# finally perform baking
				bpy.ops.object.bake_image()
				img = bpy.data.images[image.name]
				img.filepath = 'C:/Users/Camille/Documents/bcga-examples-master/imagename.png'
				img.file_format = 'PNG'
				img.save()
				# delete the high poly object and its mesh
				context.view_layer.objects.active = highPolyObject
				mesh = highPolyObject.data
				bpy.ops.object.delete()
				bpy.data.meshes.remove(mesh)
				lowPolyObject=bpy.data.objects['BCGA.002'] #added by cam 
				context.view_layer.objects.active = lowPolyObject
				# assign the baked texture to the low poly object
				blenderTexture = bpy.data.textures.new(lowPolyObject.name, type = "IMAGE")
				blenderTexture.image = image
				blenderTexture.use_alpha = True
				mat = bpy.data.materials.new(lowPolyObject.name)
				#textureSlot = material.texture_slots.add()
				#textureSlot.texture = blenderTexture
				#textureSlot.texture_coords = "UV"
				#textureSlot.uv_layer = "bcga"
				mat.use_nodes = True
				bsdf = mat.node_tree.nodes["Principled BSDF"]
				texImage = mat.node_tree.nodes.new('ShaderNodeTexImage')
				texImage.image = image
				mat.node_tree.links.new(bsdf.inputs['Base Color'], texImage.outputs['Color'])
				lowPolyObject.data.materials.append(mat)
		return {"FINISHED"}


#### Create operator that creates footprints #### 
class FootprintSet(bpy.types.Operator):
	bl_idname = "object.footprint_set"
	bl_label = "BCGA footprint"
	bl_description = "Set a building footprint for BCGA"
	bl_options = {"REGISTER", "UNDO"}
	
	# create custom property with the different sizes for the footprint 
	size: bpy.props.EnumProperty(
		items = [
			("35x15", "rectangle 35x15", "35x15"),
			("20x10", "rectangle 20x10", "20x10"),
			("10x10", "rectangle 10x10", "10x10")
		]
	)
	
	def execute(self, context):
		lightOffset = 20
		lightHeight = 20
		scene = context.scene
		# delete active object if it is a mesh
		active = context.object
		if active and active.type=="MESH":
			bpy.ops.object.delete()
		# getting width and height of the footprint
		w, h = [float(i) for i in self.size.split("x")]
		# add lights
		rx = math.atan((h+lightOffset)/lightHeight)
		rz = math.atan((w+lightOffset)/(h+lightOffset))
		def lamp_add(x, y, rx, rz):
			bpy.ops.object.light_add(
				type="SUN",
				location=((x,y,lightHeight)),
				rotation=(rx, 0, rz)
			)
			context.active_object.data.energy = 0.5
		lamp_add(w+lightOffset, h+lightOffset, -rx, -rz)
		lamp_add(-w-lightOffset, h+lightOffset, -rx, rz)
		lamp_add(-w-lightOffset, -h-lightOffset, rx, -rz)
		lamp_add(w+lightOffset, -h-lightOffset, rx, rz)
		
		# create the footprint rectangle 
		create_rectangle(context, w, h)
		align_view(context.object)
		
		return {"FINISHED"}


class FirstEdgeYmin(bpy.types.Operator):
	bl_idname = "object.first_edge_ymin"
	bl_label = "Contains min Y"
	bl_description = "The first edge contains the vertex with minimal Y coordinate and has a longer length"
	bl_options = {"REGISTER", "UNDO"}
	
	def execute(self, context):
		first_edge_ymin(context)
		return {"FINISHED"}


classes = (
	CustomColorProperty,
	CustomFloatProperty,
    FirstEdgeYmin,
    FootprintSet,
    Bake,
	Pro,
	ProMainPanel,
	BakingPanel,
	FirstEdgePanel
)
register, unregister = bpy.utils.register_classes_factory(classes)
