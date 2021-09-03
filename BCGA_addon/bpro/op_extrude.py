import pro
from pro import context
from .op_decompose import decompose_execute

class Extrude(pro.op_extrude.Extrude):
	def execute(self):
		# get original 2D shape (get State() in pro.__init__)
		state = context.getState()
		# extrude the 2D footprint using the extrude fction from pro.shape.py, returns a 3D shape 
		shape = state.shape.extrude(self)

		# if we specified the different parts of the extruded shape into the extrude fct, calls decompose_execute in bpro.op_decompose_execute
		if self.parts:
			decompose_execute(shape, self.parts)
		else:
			state.shape = shape
	
	def execute_join(self, band):
		band.extrude()