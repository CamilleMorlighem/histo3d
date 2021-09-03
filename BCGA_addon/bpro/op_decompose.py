import pro
from pro import context

def decompose_execute(shape, _parts):	
	# create a dict from the operatorDef list
	parts = {}
	for part in _parts:
		parts[part.value] = part
	# calls decompose fct in shape.py which splits the 3D shape into different labelled surface based on the surface normal 
	# components is a dictionary with a comp-selector as the key and a list of 2D-shapes as the related value 
	components = shape.decompose(parts)
	if len(components)>0:
		# now apply the rule for each decomposed 2D-shape (i think it executes each rule specified in each part ==> go on with the other rules in the rulefile hierarchically)
		for selector in components:
			for _shape in components[selector]:
				context.pushState(shape=_shape)
				parts[selector].execute()
				context.popState()


class Decompose(pro.op_decompose.Decompose):
	
	def execute(self):
		decompose_execute(context.getState().shape, self.parts)