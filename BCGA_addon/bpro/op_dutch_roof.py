import pro
import math 
from pro import context
from .polygon import Roof
from .polygon_manager import Manager
from pro import base
import math 
import numpy as np
import mathutils


class DutchRoof(pro.op_dutch_roof.DutchRoof):
    def execute(self):
        shape = context.getState().shape
        face = shape.face
        self.init(len(face.verts))
        manager = Manager()
        roof = Roof(face.verts, shape.getNormal(), manager)
        # roofheights 
        context.bm.faces.ensure_lookup_table()
        # z coordinate of the roof top is the ground building elevation + the height value 
        top_height_z = context.bm.faces[0].verts[0].co.z + self.top_height 
        # roof_height is difference between z value of the top and z value of the base of roof 
        self.roof_height= top_height_z - face.verts[0].co.z
       
        # soffits
        if self.soffits:
            manager.rule = self.soffit
            roof.inset(*self.soffits, negate=True)
        # fascias
        if self.fasciaSize:
            manager.rule = self.fascia
            roof.translate(self.fasciaSize)
        # step height and nber of steps 
        #if user provided n_steps, compute step_height 
        if self.step_height == None: 
            self.step_height= self.roof_height/self.n_steps.value
    
        #if step_height was provided, recompute it so that roof_height % step_height ==0
        else: 
            n_steps=self.roof_height//self.step_height.value
            self.step_height=self.roof_height/n_steps 

        self.roof_height=self.roof_height - (2*self.step_height)

        #find edge index based on mid_pt edge given by user 
        e_idx = 0
        for e in face.edges: 
            half_x = (e.verts[0].co.x + e.verts[1].co.x)/2
            half_y = (e.verts[0].co.y + e.verts[1].co.y)/2
            #half_z = (e.verts[0].co.z + e.verts[1].co.z)/2
            if ((half_x - self.edge_midpt[0] < 0.0001) and (half_x - self.edge_midpt[0] > -0.0001)) and ((half_y- self.edge_midpt[1] < 0.0001) and (half_y - self.edge_midpt[1] > -0.0001)): 
                break 
            e_idx = e_idx + 1 
       
        # pitches 
        pitches = []
        for i in range(len(face.edges)): 
            pitches.append(i)
        nedge = len(pitches)

        #find ajdacent edges of the facade edge 
        if e_idx >= len(face.edges): 
            e_idx = 0
        current_edge = face.edges[e_idx]
        prev_edge = face.edges[e_idx-1]
        if e_idx == len(face.edges)-1: 
            next_edge = face.edges[0]
        else: 
            next_edge = face.edges[e_idx+1]
        
        #compute intersection betw the ajdacent roof face and the roof ground face so that the z coord of the intersection point is at top height z 
        #pitch depends on roof height and "width" of the footprint 
        if prev_edge.verts[0] != current_edge.verts[0] and prev_edge.verts[0] != current_edge.verts[1]: 
            x1, y1, z1 = prev_edge.verts[0].co.x, prev_edge.verts[0].co.y, prev_edge.verts[0].co.z 
            x2, y2, z2 = prev_edge.verts[1].co.x, prev_edge.verts[1].co.y, prev_edge.verts[1].co.z 
        else:
            x1, y1, z1 = prev_edge.verts[1].co.x, prev_edge.verts[1].co.y, prev_edge.verts[1].co.z 
            x2, y2, z2 = prev_edge.verts[0].co.x, prev_edge.verts[0].co.y, prev_edge.verts[0].co.z
        
        x3, y3, z3 = next_edge.verts[1].co.x, next_edge.verts[1].co.y, next_edge.verts[1].co.z 
        x4, y4, z4 = next_edge.verts[0].co.x, next_edge.verts[0].co.y, next_edge.verts[0].co.z 

        a= np.array([[-(y1-y2), -(x3-x4)], [(x1-x2), -(y3-y4)]])
        #a= np.array([[(x4-x2), -(x3-x4)], [(y4-y2), -(y3-y4)]])
        b = np.array([x4-x2, y4-y2])
        coeff = np.linalg.solve(a,b)
        parameters=list(coeff)
        intersect_pt = mathutils.Vector((parameters[1]*(x3-x4) + x4, parameters[1]*(y3-y4) + y4, z4))
        base_length = math.sqrt((x2 - intersect_pt.x)**2 + (y2 - intersect_pt.y)**2 +(z2 - intersect_pt.z)**2)
        
        #pitches of the long edges are set to 90° (the stepped facade is not pitched in old dutch roofs)
        #pitches of the short edges are computed using the roof_height and edge length to make sure top of the roof is flat ! 
        #use edge index to order the pitches values 
        #edge index corresponds to facade edge ==> pitches computed using roof heigh and edge length, while the other adjacent edges have 90° pitches
        a=0
        for i in range(len(face.edges)): 
            if i == len(face.edges)-1:
                next_i=0
            else : 
                next_i=i+1

            if i == 0: 
                pitches[e_idx-nedge+i] = 90
            elif i == 2: 
                pitches[e_idx-nedge+i] = self.pitches[0]
                #pitches[e_idx-nedge+i] = 90
            else: 
                if a == 0: 
                    j= math.degrees(math.atan(self.roof_height/(base_length/2))) 
                    pitches[e_idx-nedge+i] = j #math.degrees(math.atan(self.roof_height/(face.edges[e_idx-nedge+next_i].calc_length()/2))) 
                #pitches[e_idx-nedge+i] = 30
                    a=a+1
                else: 
                    pitches[e_idx-nedge+i] = j 

        self.pitches=pitches
        
        """
        pitches = []
        j=0
        for i in range(len(face.edges)): 
            if i == len(face.edges)-1:
                next_i=0
            else : 
                next_i=i+1

            if face.edges[i].calc_length() < face.edges[next_i].calc_length(): 
                #pitch=math.degrees(math.atan(self.roof_height/(face.edges[i].calc_length()/2))) 
                #pitches.append(pitch)
                if j==0: 
                    pitches.append(90)
                    j=j+1
                else: 
                    pitches.append(self.pitches[0])
            else: 
                #pitches.append(90)
                pitch=math.degrees(math.atan(self.roof_height/(face.edges[next_i].calc_length()/2))) 
                pitches.append(pitch)

        
        self.pitches=pitches
        """
        # hip roof itself
        manager.rule = self.face
        roof.stepped_roof(shape, top_height_z, self.step_height, self.pitches, self.facade_thickness, e_idx)
        shape.delete()
        # finalizing: if there is a rule for the shape, execute it
        for entry in manager.shapes:
            context.pushState(shape=entry[0])
            entry[1].execute()
            context.popState()
        
        