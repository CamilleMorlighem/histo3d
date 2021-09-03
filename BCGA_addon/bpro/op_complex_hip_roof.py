import pro
from pro import context
from .polygon import Roof
from .polygon_manager import Manager

class ComplexHipRoof(pro.op_complex_hip_roof.ComplexHipRoof):
    def execute(self):
        shape = context.getState().shape
        face = shape.face
        self.init(len(face.verts))
        manager = Manager()
        roof = Roof(face.verts, shape.getNormal(), manager)
        # roofheight 
        context.bm.faces.ensure_lookup_table()
        # z coordinate of the roof top is the ground building elevation + the height value 
        top_height_z = context.bm.faces[0].verts[0].co.z + self.top_height 
       
        # soffits
        if self.soffits:
            manager.rule = self.soffit
            roof.inset(*self.soffits, negate=True)
        # fascias
        if self.fasciaSize:
            manager.rule = self.fascia
            roof.translate(self.fasciaSize)
        # hip roof itself
        manager.rule = self.face
        roof.complex_roof(shape, top_height_z)
        #if you comment the following line, the top horizontal face will be part of the output object (if u uncomment it, then you have no separation between the building and the roof )
        shape.delete()
        # finalizing: if there is a rule for the shape, execute it
        for entry in manager.shapes:
            context.pushState(shape=entry[0])
            entry[1].execute()
            context.popState()
        