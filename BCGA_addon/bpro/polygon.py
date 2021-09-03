import math
import bpy, bmesh
from pro import context
import mathutils 
import numpy as np 
import copy 

from .shape import createShape2d, createRectangle, Shape2d
from .internal import normalOfPolygon
from .internal import nearestPointOfLines, straightSkeletonOfPolygon
from .util import zero
from bpro.op_extrude import Extrude 



timeTolerance = 0.0001
distanceTolerance2 = 0.000001

class Manager:
    """A dummy manager for Polygon class declared below"""
    
    def getValue(self, obj):
        # return obj as is
        return obj
    
    def resolve(self, shape, value):
        # do nothing here
        pass
    
# create an instance of the dummy manager class
_manager = Manager()
    

class Polygon:
    def __init__(self, verts, axis, manager=None):
        self.axis = axis
        self.manager = manager if manager else _manager
        edges = []
        corners = []
        numVerts = len(verts)
        self.numEdges = numVerts
        numVerts -= 1
        i = 0
        vert = verts[0]
        vec = vert.co - verts[numVerts].co
        vec.normalize()
        _edge = Edge(vec, axis)
        prevEdge = _edge
        while i < numVerts:
            vert2 = verts[i+1]
            vec = vert2.co-vert.co
            vec.normalize()
            edge = Edge(vec, axis)
            edges.append(edge)
            corners.append(Corner(vert, prevEdge, edge, axis))
            vert = vert2
            prevEdge = edge
            i += 1
        edges.append(_edge)
        corners.append(Corner(vert, prevEdge, _edge, axis))
        self.edges = edges
        self.corners = corners
        # assign a pair of edges to each corner
        nextCorner = corners[0]
        # i = numVerts from the previous while cycle
        while i>=0:
            edge = edges[i]
            corner = corners[i]
            edge.corner1 = corner
            edge.corner2 = nextCorner
            nextCorner = corner
            i -= 1
    
    def inset(self, *distances, **kwargs):
        manager = self.manager
        translate = kwargs["height"]*self.axis if "height" in kwargs else None
        negate = kwargs["negate"] if "negate" in kwargs else False
        corners = self.corners
        distancePerEdge = False if len(distances)==1 else True
        #if there is only 3 inset values, the 2 last corners have the same insets 
        if distancePerEdge:
            distance1 = manager.getValue(distances[-2])
            _d = distance1
            distance2 = manager.getValue(distances[-1])
        else:
            distance1 = manager.getValue(distances[0])
            distance2 = distance1
        corner = corners[-1]
        prevVert1 = corner._vert
        _vert1 = prevVert1
        corner.inset(distance1, distance2, translate, negate)
        prevVert2 = corner._vert
        _vert2 = prevVert2
        i = 0
        numCorners = self.numEdges-1
        while i<numCorners:
            corner = corners[i]
            if distancePerEdge:
                distance1 = distance2
                distance2 = manager.getValue(distances[i])
            vert1 = corner._vert
            corner.inset(distance1, distance2, translate, negate)
            vert2 = corner._vert
            if distance1!=0:
                manager.resolve(
                    createShape2d((prevVert1, vert1, vert2, prevVert2)),
                    distances[i-1] if distancePerEdge else distances[0]
                )
            prevVert1 = vert1
            prevVert2 = vert2
            i += 1
        if not distancePerEdge or _d!=0:
            manager.resolve(
                createShape2d((prevVert1, _vert1, _vert2, prevVert2)),
                distances[-2] if distancePerEdge else distances[0]
            )
    
    def getShape(self, constructor):
        """Returns a polygon composed of the present vertices"""
        verts = []
        for corner in self.corners:
            verts.append(corner._vert)
        face = context.bm.faces.new(verts)
        return constructor(face.loops[0])
            
    def straightSkeleton(self, getVert=None):
        sequences = {}
        seq = Sequence(self.edges[0], len(self.corners), self.axis, getVert)
        sequences[seq.id] = seq
        numSequences = 1
        
        while numSequences:
            changes = {}
            for _id in sequences:
                seq = sequences[_id]
                _sequences = seq.process(changes)
                changes[seq.id] = None
                    
            for _id in changes:
                if changes[_id]:
                    sequences[_id] = changes[_id]
                    numSequences += 1
                else:
                    del sequences[_id]
                    numSequences -= 1

    def translate(self, distance, axis=None):
        translate = distance*(axis if axis else self.axis)
        corners = self.corners
        corner = corners[-1]
        prevVert1 = corner._vert
        _vert1 = prevVert1
        corner.vert = corner.vert + translate
        corner._vert = context.bm.verts.new(corner.vert)
        prevVert2 = corner._vert
        _vert2 = prevVert2
        i = 0
        numCorners = self.numEdges-1
        while i<numCorners:
            corner = corners[i]
            vert1 = corner._vert
            corner.vert = corner.vert + translate
            corner._vert = context.bm.verts.new(corner.vert)
            vert2 = corner._vert
            self.manager.resolve(createRectangle((prevVert1, vert1, vert2, prevVert2)))
            prevVert1 = vert1
            prevVert2 = vert2
            i += 1
        self.manager.resolve(createRectangle((prevVert1, _vert1, _vert2, prevVert2)))
    
    def create_stepped_facade(self, starting_triang_face, n_steps, top_height_z):
        stepped_face = [] 
        #STEP2: find the 2 base points of that face and the top vertex 
        top_pt=starting_triang_face[0]
        base_pts=[]
        base_pts_copy=[]
        for pt in starting_triang_face: 
            if pt.co.z > top_pt.co.z: 
                top_pt=pt 
            else: 
                base_pts.append(pt)
                base_pts_copy.append(pt)
                stepped_face.append(pt)
        
        #for the whole procedure we don't use the top_pt of the triangle but we offset it by setting it at a higher height 
        #this is done to avoid self intersection in the output stepped face because of coordinates rounding 
        top_pt_offset = mathutils.Vector((top_pt.co.x, top_pt.co.y, top_pt.co.z + 0.2))

        #recompute roof_height and step_height 
        roof_height = top_pt_offset.z - base_pts[0].co.z  
        #n_steps=roof_height//step_height
        if n_steps < 1: 
            n_steps = 1
        step_height=roof_height/n_steps
      

        #create steps along 2 sides of the triangular face until we reach the top vertex of the triangle 
        while top_pt_offset.z - (base_pts[0].co.z + step_height) > 0.002: 
            #STEP3: get a', b', the base points offset by step_height 
            offset_base_pts=[]
            
            #add a',b' to the mesh vertices 
            for base in base_pts: 
                offset_base=context.bm.verts.new((base.co.x, base.co.y, base.co.z+ step_height)) 
                offset_base_pts.append(offset_base)
                stepped_face.append(offset_base)
                
            #STEP4: intersect line segment a'b' with segments ac and bc, the two sides of the triangle passing by the top vertex  
            dir_offsets=offset_base_pts[0].co - offset_base_pts[1].co

            #get param equation of the line segments and find their point of intersection 
            new_base_pts=[]
            for i in range(0, len(base_pts)): 
                dir_base_top= top_pt_offset - base_pts[i].co 
                a = np.array([[dir_base_top.x, -dir_offsets.x], [dir_base_top.y, -dir_offsets.y], [dir_base_top.z, -dir_offsets.z]])
                b = np.array([offset_base_pts[0].co.x - top_pt_offset.x, offset_base_pts[0].co.y - top_pt_offset.y, offset_base_pts[0].co.z - top_pt_offset.z])
                #use least square because we have 3 equations for two unknowns 
                param= np.linalg.lstsq(a, b, rcond=None)
                parameters=list(param[0])
                x=top_pt_offset.x + parameters[0]*dir_base_top.x
                y=top_pt_offset.y + parameters[0]*dir_base_top.y
                z=top_pt_offset.z + parameters[0]*dir_base_top.z
               
                v=context.bm.verts.new((x,y, z)) 
                stepped_face.append(v)
                new_base_pts.append(v)
            
            #repeat the code in the while loop with the intersection points now set as base points 
            base_pts=new_base_pts

        #when we are at the top of the triangular face, we create the final part at the top of the steps 
        #we add 3 last vertices (*2), the last step two vertices and one vertex of the top part 
        for base in base_pts: 
            offset_base=context.bm.verts.new((base.co.x, base.co.y, base.co.z+ step_height)) 
            half_pt_x= (offset_base.co.x + top_pt_offset.x)/2
            half_pt_y= (offset_base.co.y + top_pt_offset.y)/2
            half_pt_z= offset_base.co.z
            
            half_pt = context.bm.verts.new((half_pt_x, half_pt_y, half_pt_z))
            half_pt_offset= context.bm.verts.new((half_pt_x, half_pt_y, top_height_z)) #half_pt_z + step_height))
            stepped_face.append(offset_base)
            stepped_face.append(half_pt)
            stepped_face.append(half_pt_offset)

        #all vertices of the face are now stored into stepped_face, we need to reorder them 
        stepped_face_part1=[]
        stepped_face_part2=[]
        for i in range(0, len(stepped_face)-6): 
            if i%2==0: 
                stepped_face_part1.append(stepped_face[i])
            else: 
                stepped_face_part2=[stepped_face[i]]+stepped_face_part2
        
        #we create two faces, one face is the stepped face visible part + the triangular face of the roof 
        #the other face is only the visible part 
        ordered_face=stepped_face_part1 + stepped_face[-6:-3] + [stepped_face[-1]] + [stepped_face[-2]] + [stepped_face[-3]] + stepped_face_part2
        ordered_face_ext = [top_pt] + stepped_face_part1 + stepped_face[-6:-3] + [stepped_face[-1]] + [stepped_face[-2]] + [stepped_face[-3]] + stepped_face_part2
       
        #check faces are properly oriented (ccw)
        vtex_lst=[]
        for vtex in starting_triang_face: 
            vtex_lst.append(vtex.co)
        
        ordered_face_vtx = []
        for vtex in ordered_face: 
            ordered_face_vtx.append(vtex.co)

        #if the starting triangular face and full stepped face have opposite normals, the stepped face is badly oriented 
        if normalOfPolygon(ordered_face_vtx).dot(normalOfPolygon(vtex_lst)) < 0: 
            ordered_face.reverse()
        #the face containing the visible part must be oriented in the opposite direction than the full stepped face because there are on opposite sides of the mesh 
        else: 
            ordered_face_ext.reverse()
        
        return (ordered_face, ordered_face_ext)
        


class Roof(Polygon):
    def __init__(self, verts, axis, manager=None):
        super().__init__(verts, axis, manager)

    def roof(self, *pitches):
        manager = self.manager 
        numPitches = len(pitches)
        if numPitches==1:
            pitch = pitches[0]
            distance = 1/math.tan( math.radians(manager.getValue(pitch)) )
            for edge in self.edges:
                edge.pitch = pitch
                edge.distance = distance
        else:
            # each edge has its own pitch
            i = 0
            while i<numPitches:
                edge = self.edges[i]
                pitch = pitches[i]
                edge.pitch = pitch
                # edge.distance is equal to cotangent of pitch
                edge.distance = 1/math.tan( math.radians(manager.getValue(pitch)) )
               
                i += 1
        def getVert(vert, t):
            return vert + t*self.axis
        self.straightSkeleton(getVert)

        for edge in self.edges:
            edge.leftVerts.reverse()
            face = [edge.leftVerts.pop()]
            face += edge.rightVerts
            face += edge.leftVerts
            
            shape = createShape2d(face)
            self.manager.resolve(shape, edge.pitch)


    def gable_hip_roof(self, base_shape, top_height_z, pitches, edge_index): 
        manager = self.manager
        face=base_shape.face
        #create list with vertices of base roof shape 
        lst_vertex=[]
        for v in face.verts: 
            lst_vertex.append(v.co)
        
        numPitches = len(pitches)

        # each edge has its own pitch
        i = 0
        while i<numPitches:
            edge = self.edges[i]
            pitch = pitches[i]
            edge.pitch = pitch
            # edge.distance is equal to cotangent of pitch
            edge.distance = 1/math.tan( math.radians(manager.getValue(pitch)) )
            i += 1
        def getVert(vert, t):
            return vert + t*self.axis
        self.straightSkeleton(getVert)

        starting_triang_face = None 
        top_pts_list = set()
        edge_i = 0
        for edge in self.edges:
            edge.leftVerts.reverse()
            face = [edge.leftVerts.pop()]
            face += edge.rightVerts
            face += edge.leftVerts
            shape = createShape2d(face)
            self.manager.resolve(shape, edge.pitch)
            #for z coordinate of the roof top vertices 
            if edge_i == edge_index : 
                starting_triang_face = face 
            for pt in face: 
                if pt.co not in lst_vertex: 
                    top_pts_list.add((pt.co.z, pt.co.y, pt.co.x))
            edge_i += 1

        coords = False 
        coor = False 
        
        #following needed to ensure the proper roof height (needed for slanted line roof top)
        #if the face which matches with edge index has 3 vertices
        if len(self.edges) == 4 and starting_triang_face != None and len(starting_triang_face)==3: 
            #get top pt and base pts 
            face = starting_triang_face
            top_pt=face[0]
            base_pts_copy=[]
            for pt in face: 
                if pt.co.z > top_pt.co.z: 
                    top_pt=pt 
                else: 
                    base_pts_copy.append(pt)

            #find the top vertices of the roof 
            for edge in top_pt.link_edges:
                if (edge.verts[0].co != base_pts_copy[0].co)  and (edge.verts[0].co != base_pts_copy[1].co) and (edge.verts[1].co != base_pts_copy[0].co) and (edge.verts[1].co != base_pts_copy[1].co): 
                    #find the highest and lowest 
                    if edge.verts[0].co.z > edge.verts[1].co.z: 
                        highest = edge.verts[0]
                        lowest = edge.verts[1]
                    else: 
                        highest = edge.verts[1]
                        lowest = edge.verts[0]
                    x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                    x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z

                    #move the highest along its adjacent trinagular face (so that preserve the pitch value) until it has the correct roofheight 
                    for fa in highest.link_faces: 
                        x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                        x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z
                        #find intersection point betw triangular face and translated line which intersects the face in z = top height z 
                        if len(fa.verts)==3: 
                            top_pt=fa.verts[0]
                            base_pts_highest=[]
                            for pt in fa.verts: 
                                if pt.co.z > top_pt.co.z: 
                                    top_pt=pt 
                                else: 
                                    base_pts_highest.append(pt)
                            fa.normal_update()
                           
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                            x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z 
                            a= np.array([[1,0, -(x1-x2), 0], [0,1, -(y1-y2), 0],  [0,0, -(z1-z2), -1], [n1, n2, 0, 0]])
                            b = np.array([x2, y2, -top_height_z, n1*x2 + n2*y2 - n3*(top_height_z - z2)])
                            
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord1 = mathutils.Vector((parameters[0], parameters[1], top_height_z))
                            coor = True 

                    #now we need to move the lowest top vertex 
                    for fa in lowest.link_faces: 
                        if len(fa.verts)==3 and coor: 
                            fa.normal_update()
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            '''
                            a = np.array([[1,0,0, -(x1-x2)], [0,1,0, -(y1-y2)], [0,0,1, -(z1-z2)], [n1, n2, n3, 0]])
                            b = np.array([x2, y2, parameters[3], n1*x1 + n2*y1 + n3*z1])
                           
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord2 = mathutils.Vector((parameters[0], parameters[1], parameters[2]))
                            coords = True 
                            '''
                            #find intersection betw opposite triangular face and the 2 other roof faces that all pass through the previous computed point and base points of the roof base 
                            planes = []
                            for x in base_pts_highest: 
                                for e in x.link_edges:
                                    if e.verts[0] not in base_pts_highest or e.verts[1] not in base_pts_highest:  
                                        if e.verts[0].co != x.co and e.verts[0].co.z == x.co.z: 
                                            v1, v2, v3 = normalOfPolygon([x.co, e.verts[0].co, new_coord1]).normalized()
                                            #plane.append([x.co, e.verts[0].co, new_coord1])
                                            gp1 = e.verts[0].co
                                            planes.append((gp1, (v1, v2, v3)))
                                        elif e.verts[1].co !=x.co and e.verts[1].co.z == x.co.z: 
                                            #plane.append([x.co, e.verts[1].co, new_coord1])
                                            v1,v2,v3 = normalOfPolygon([x.co, e.verts[1].co, new_coord1]).normalized()
                                            gp1 = e.verts[1].co
                                            planes.append((gp1, (v1, v2, v3)))
                                            
                                        
                            v1, v2, v3 = planes[0][1]
                            gp1 = planes[0][0]
                            r1, r2, r3 = planes[1][1]
                            gp2 = planes[1][0]
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            
                            a = np.array([[n1, n2, n3], [v1, v2, v3], [r1, r2, r3]])
                            b = np.array([n1*gp1.x + n2*gp1.y + n3*gp1.z, v1*gp1.x + v2*gp1.y + v3*gp1.z, r1*gp2.x + r2*gp2.y + r3*gp2.z])
                            
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord2 = mathutils.Vector((parameters[0], parameters[1], parameters[2]))
                            coords = True 

                    #replace coords    
                    if coords  and coor: 
                        highest.co = new_coord1
                        lowest.co = new_coord2
                        
       


    def stepped_roof(self, base_shape, top_height_z, step_height, pitches, facade_thickness, edge_index):
        if facade_thickness==None: 
            facade_thickness= 0.5

        manager = self.manager
        face=base_shape.face
        #create list with vertices of base roof shape 
        lst_vertex=[]
        for v in face.verts: 
            lst_vertex.append(v.co)

        numPitches = len(pitches)
        roof_height_z = top_height_z - (2*step_height)

        #STEP0: make basic hip roof 
        # each edge has its own pitch
        i = 0
        while i<numPitches:
            edge = self.edges[i]
            pitch = pitches[i]
            edge.pitch = pitch
            # edge.distance is equal to cotangent of pitch
            edge.distance = 1/math.tan( math.radians(manager.getValue(pitch)) )
            i += 1
        def getVert(vert, t):
            return vert + t*self.axis
        self.straightSkeleton(getVert)

        highest_z = 0 
        #STEP1: find one triangular roof face where to built the stepped roof (which has index edge_idx)
        starting_triang_face=None
        edge_i = 0
        top_pts_list = set()
        faces = []
        for edge in self.edges:
            edge.leftVerts.reverse()
            face = [edge.leftVerts.pop()]
            face += edge.rightVerts
            face += edge.leftVerts   
            if edge_i == edge_index: #len(face)==3 and starting_triang_face==None: #edge_i == edge_index: 
                starting_triang_face=face 
            #for all the other faces of the straight skeleton, add them to the mesh 
            else: 
                shape = createShape2d(face)
                self.manager.resolve(shape, edge.pitch)
            edge_i=edge_i+1
            faces.append(face)
            #keep track of points at the roof top to later force the roofheight 
            for pt in face: 
                if pt.co not in lst_vertex: 
                    top_pts_list.add((pt.co.x, pt.co.y, pt.co.z))
                    if highest_z < pt.co.z: 
                        highest = pt.co 
                        highest_z = pt.co.z 
               
        #find base points and top point of triangle 
        top_pt=starting_triang_face[0]
        base_pts_copy=[]
        for pt in starting_triang_face: 
            if pt.co.z > top_pt.co.z: 
                top_pt=pt 
            else: 
                base_pts_copy.append(pt)
        
        shape = createShape2d(starting_triang_face)
        self.manager.resolve(shape)
        shape.face.normal_update()
        #### 
        coords = False 
        coor = False
       
        if (highest_z - top_pt.co.z) > 0.2 : 
            new_height = top_height_z 
        else: 
            new_height = top_height_z - 2 *step_height
        #following needed to ensure the proper roof height (needed for slanted line roof top)
        #if the face which matches with edge index has 3 vertices
        if len(self.edges) == 4 and starting_triang_face != None and len(starting_triang_face)==3: 
            #get top pt and base pts 
            face = starting_triang_face
            top_pt=face[0]
            base_pts_copy=[]
            for pt in face: 
                if pt.co.z > top_pt.co.z: 
                    top_pt=pt 
                else: 
                    base_pts_copy.append(pt)

            #find the top vertices of the roof 
            for edge in top_pt.link_edges:
                if (edge.verts[0].co != base_pts_copy[0].co)  and (edge.verts[0].co != base_pts_copy[1].co) and (edge.verts[1].co != base_pts_copy[0].co) and (edge.verts[1].co != base_pts_copy[1].co): 
                    #find the highest and lowest 
                    if edge.verts[0].co.z > edge.verts[1].co.z: 
                        highest = edge.verts[0]
                        lowest = edge.verts[1]
                    else: 
                        highest = edge.verts[1]
                        lowest = edge.verts[0]
                    x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                    x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z

                    #move the highest along its adjacent trinagular face (so that preserve the pitch value) until it has the correct roofheight 
                    for fa in highest.link_faces: 
                        x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                        x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z
                        #find intersection point betw triangular face and translated line which intersects the face in z = top height z 
                        if len(fa.verts)==3: 
                            top_pt=fa.verts[0]
                            base_pts_highest=[]
                            for pt in fa.verts: 
                                if pt.co.z > top_pt.co.z: 
                                    top_pt=pt 
                                else: 
                                    base_pts_highest.append(pt)
                            fa.normal_update()
                        
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            x1, y1, z1 = lowest.co.x, lowest.co.y, lowest.co.z 
                            x2, y2, z2 = highest.co.x, highest.co.y, highest.co.z 
                            a= np.array([[1,0, -(x1-x2), 0], [0,1, -(y1-y2), 0],  [0,0, -(z1-z2), -1], [n1, n2, 0, 0]])
                            b = np.array([x2, y2, -new_height, n1*x2 + n2*y2 - n3*(new_height - z2)])
                            
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord1 = mathutils.Vector((parameters[0], parameters[1], new_height))
                            coor = True 

                    #now we need to move the lowest top vertex 
                    for fa in lowest.link_faces: 
                        if len(fa.verts)==3 and coor: 
                            fa.normal_update()
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            
                            #find intersection betw opposite triangular face and the 2 other roof faces that all pass through the previous computed point and base points of the roof base 
                            planes = []
                            for x in base_pts_highest: 
                                for e in x.link_edges:
                                    if e.verts[0] not in base_pts_highest or e.verts[1] not in base_pts_highest:  
                                        if e.verts[0].co != x.co and e.verts[0].co.z == x.co.z: 
                                            v1, v2, v3 = normalOfPolygon([x.co, e.verts[0].co, new_coord1]).normalized()
                                            #plane.append([x.co, e.verts[0].co, new_coord1])
                                            gp1 = e.verts[0].co
                                            planes.append((gp1, (v1, v2, v3)))
                                        elif e.verts[1].co !=x.co and e.verts[1].co.z == x.co.z: 
                                            #plane.append([x.co, e.verts[1].co, new_coord1])
                                            v1,v2,v3 = normalOfPolygon([x.co, e.verts[1].co, new_coord1]).normalized()
                                            gp1 = e.verts[1].co
                                            planes.append((gp1, (v1, v2, v3)))
                                            
                                        
                            v1, v2, v3 = planes[0][1]
                            gp1 = planes[0][0]
                            r1, r2, r3 = planes[1][1]
                            gp2 = planes[1][0]
                            n1, n2, n3 = fa.normal.x, fa.normal.y, fa.normal.z 
                            
                            a = np.array([[n1, n2, n3], [v1, v2, v3], [r1, r2, r3]])
                            b = np.array([n1*gp1.x + n2*gp1.y + n3*gp1.z, v1*gp1.x + v2*gp1.y + v3*gp1.z, r1*gp2.x + r2*gp2.y + r3*gp2.z])
                            
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord2 = mathutils.Vector((parameters[0], parameters[1], parameters[2]))
                            coords = True 

                    #replace coords    
                    if coords  and coor: 
                        highest.co = new_coord1
                        lowest.co = new_coord2
                    
        top_pt=starting_triang_face[0]
        base_pts_copy=[]
        for pt in starting_triang_face: 
            if pt.co.z > top_pt.co.z: 
                top_pt=pt 
            else: 
                base_pts_copy.append(pt)

        ####
        top_height_z = 2 * step_height + top_pt.co.z
        roof_height_z = top_pt.co.z
        #recompute roof_height and n steps 
        roof_height = roof_height_z - base_pts_copy[0].co.z
       
        n_steps=roof_height//step_height
      

        #STEP2: get the stepped facade of the non translated triangular face 
        facade = []
        order_face, order_face_ext = self.create_stepped_facade(starting_triang_face, n_steps, top_height_z)
        for pt in order_face:
            facade.append(mathutils.Vector((pt.co.x, pt.co.y, pt.co.z)))
        
        #store the vertices into another list before removing them from the mesh 
        for x in order_face:
            if x not in starting_triang_face: 
                bmesh.ops.delete(context.bm, geom = [x], context = 'VERTS' )

        #shape = createShape2d(starting_triang_face)
        #self.manager.resolve(shape)
        #shape.face.normal_update()
        
        #STEP3: offset the footprint towards the interior of the polygon
        #1 get normal of starting triangle 
        starting_plane_normal = mathutils.Vector((shape.face.normal.x, shape.face.normal.y, shape.face.normal.z))
        #2 get one point in the plane parallel to the starting triangle at a distance of facade thickness twoards the interior of the fp 
        point_in_plane = False 
        for point in shape.face.verts: 
            for e in point.link_edges: 
                if not (e.verts[0] in shape.face.verts and e.verts[1] in shape.face.verts) and len(point.link_edges)==4: 
                    #get point belonging to the footprint 
                    if e.verts[0].co.z != e.verts[1].co.z: 
                        if e.verts[0] != point: 
                            point_in_plane = e.verts[0]
                        else: 
                            point_in_plane = e.verts[1]
                        ground_height = point_in_plane.co.z

            if point_in_plane: 
                #compute coordinate of point in the parallel plane at a distance equals to facade_thickness 
                point_in_plane = point_in_plane.co + (-1) * starting_plane_normal * facade_thickness 
                break 

        #compute the intersection between the roof edges and the parallel plane 
        new_pt_coords = []
        pts_to_be_moved = []
        
        n1, n2, n3 = starting_plane_normal.x, starting_plane_normal.y, starting_plane_normal.z 
        p1, p2, p3 = point_in_plane.x, point_in_plane.y, point_in_plane.z 
        #for each point in the triangular face 
        for point in shape.face.verts: 
            direction = False 
            #for each edge ajdacent to the point 
            for e in point.link_edges: 
                #if the edge is not a side of the triangular face 
                if not (e.verts[0] in shape.face.verts and e.verts[1] in shape.face.verts):
                    #if there are 4 adjacent edges (we are on the base points of the triangle)
                    if len(point.link_edges)==4: 
                        #if the z coordinates of the edge are the same, the edge is horizontal (roof edge)
                        if e.verts[0].co.z == e.verts[1].co.z: 
                            pts_to_be_moved.append(point) 
                            #check intersection between line equation of edge and plane 
                            x1, y1, z1 = e.verts[0].co.x, e.verts[0].co.y, e.verts[0].co.z 
                            x2, y2, z2 = e.verts[1].co.x, e.verts[1].co.y, e.verts[1].co.z 

                            a = np.array([[1,0,0, -(x1-x2)], [0,1,0, -(y1-y2)], [0,0,1, -(z1-z2)], [n1, n2, n3, 0]])
                            b = np.array([x2, y2, z2, n1*p1 + n2*p2 + n3*p3])
                            coeff = np.linalg.solve(a,b)
                            parameters=list(coeff)
                            new_coord = mathutils.Vector((parameters[0], parameters[1], parameters[2]))
                           
                            new_pt_coords = new_pt_coords + [new_coord]
                           
                            if direction == False: 
                                new_pt_coords = new_pt_coords + [mathutils.Vector((new_coord.x, new_coord.y, ground_height))]
                                direction = True 
                            else: 
                                new_pt_coords[-2] = mathutils.Vector((new_coord.x, new_coord.y, ground_height))
                                
                        #if the z coordiantes is not constant, edge is vertical (link ground and roof)
                        else: 
                            if e.verts[0] != point: 
                                pts_to_be_moved.append(e.verts[0])
                                facade.append(mathutils.Vector((e.verts[0].co.x,e.verts[0].co.y, e.verts[0].co.z)))
                                ground_pt = e.verts[0].co
                            else: 
                                pts_to_be_moved.append(e.verts[1])
                                facade.append(mathutils.Vector((e.verts[1].co.x,e.verts[1].co.y, e.verts[1].co.z)))
                                ground_pt = e.verts[1].co
                            if direction == False:
                                new_pt_coords.append(mathutils.Vector((ground_pt.x, ground_pt.y, ground_pt.z))) 
                                direction = True 

                    #if there are 3 ajdacent edges, we are at the top edge of the roof     
                    else: 
                        pts_to_be_moved.append(point) 
                        #check intersection between line equation of edge and plane 
                        x1, y1, z1 = e.verts[0].co.x, e.verts[0].co.y, e.verts[0].co.z 
                        x2, y2, z2 = e.verts[1].co.x, e.verts[1].co.y, e.verts[1].co.z 

                        a = np.array([[1,0,0, -(x1-x2)], [0,1,0, -(y1-y2)], [0,0,1, -(z1-z2)], [n1, n2, n3, 0]])
                        b = np.array([x2, y2, z2, n1*p1 + n2*p2 + n3*p3])
                        coeff = np.linalg.solve(a,b)
                        parameters=list(coeff)
                        new_coord = mathutils.Vector((parameters[0], parameters[1], parameters[2]))
                        new_pt_coords.append(new_coord)

        #for each point of the existing facade, replace their coordinates with point located on a parallel plane translated towards the interior of the fp 
        for i in range (0, len(pts_to_be_moved)):
            pts_to_be_moved[i].co = new_pt_coords[i]
        
        #remove trinagular face from mesh 
        shape.face.normal_update()
        bmesh.ops.delete(context.bm, geom = [shape.face], context = 'FACES_ONLY' )
        #recompute top and base points (were translated)
        top_pt=starting_triang_face[0]
        base_pts_copy=[]
        for pt in starting_triang_face: 
            if pt.co.z > top_pt.co.z: 
                top_pt=pt 
            else: 
                base_pts_copy.append(pt)

        #STEP4: get stepped facade of transalted triangular face 
        ordered_face, ordered_face_ext = self.create_stepped_facade(starting_triang_face, n_steps, top_height_z)
        #add both faces to the mesh 
        shape = createShape2d(ordered_face_ext)
        self.manager.resolve(shape)
        shape.face.normal_update()

        shape = createShape2d(ordered_face)
        self.manager.resolve(shape)
        shape.face.normal_update()
       
        #STEP5: extrude the stepped face newly created 
        #first find the wall face "below" the stepped face (together they form the building façade)
        selected_faces=[]
        edges=[]
        
        for e in shape.face.edges:
            edges.append(e)
        for e in edges:
            #find the common edge to the stepped face and the wall face on the same façade 
            if (e.verts[0].co==base_pts_copy[0].co or e.verts[0].co==base_pts_copy[1].co) and (e.verts[1].co==base_pts_copy[0].co or e.verts[1].co==base_pts_copy[1].co): 
                #get all the faces that share this edge 
                linked=e.link_faces
                for fa in linked: 
                    #the wall face is the one with a normal in the same direction as the stepped_face normal 
                    fa.normal_update()
                    if shape.face.normal.dot(fa.normal) > 0.9: 
                        selected_faces.append(fa)
        
        #the full stepped face and wall face are joined to form the building facade and then the merged face is extruded
        f = bmesh.utils.face_join(selected_faces)
        f.normal_update()
        dir_translation = mathutils.Vector((f.normal.x, f.normal.y, f.normal.z)) 
        r = bmesh.ops.extrude_face_region(context.bm, geom=[f], use_normal_flip=False)
        verts = [e for e in r['geom'] if isinstance(e, bmesh.types.BMVert)]
        TranslateDirection = dir_translation * (facade_thickness) # Extrude Strength/Length
        bmesh.ops.translate(context.bm, vec = TranslateDirection, verts=verts)
        
        #remove non needed faces (output is a solid, no "faces" in the interior of the mesh)
        bmesh.ops.delete(context.bm, geom = selected_faces, context = 'FACES_ONLY' )
        
        #STEP6: iterate over the vertices of the newly extruded face and replace their coord with stepped face coords of non translated triangular face 
        context.bm.faces.ensure_lookup_table()
        found_face = False 
        
        for fa in context.bm.faces: 
            fa.normal_update()
            #if the face fa and merged face have the same normal, we may have found extruded face 
            if (fa.normal.dot(-1*f.normal) >= 0.9) and (fa.verts[0].co != f.verts[0].co): 
                #if the face contains points with z equals to top height z we found the extruded face 
                for pt in fa.verts: 
                    if ((pt.co.z - top_height_z) >= -0.5) and ((pt.co.z - top_height_z) <= 0.5)  :
                        found_face = True  
                        break 

                #iterate over its coordinates and replace them 
                if found_face: 
                    for i in range(0, len(fa.verts)): 
                        fa.verts[i].co.x = facade[i].x
                        fa.verts[i].co.y = facade[i].y
                        #fa.verts[i].co.z = facade[i].z 
                fa.normal_update()
        
        #remove non needed faces (output is a solid, no "faces" in the interior of the mesh)
        bmesh.ops.delete(context.bm, geom = [f], context = 'FACES_ONLY' )  

           
    def complex_roof(self, base_shape, top_height_z): 
        manager = self.manager
        face=base_shape.face
        #create list with vertices of base roof shape 
        lst_vertex=[]
        for v in face.verts: 
            lst_vertex.append(v.co)

        #get straight skeleton 
        verts, edges, faces = straightSkeletonOfPolygon(polygon_vertices=lst_vertex)

        top_pts_list = set()
        #for each face, get its vertices and add them to the mesh 
        for face_range in faces: 
            v_lst=[]
            for idx in face_range: 
                pt=verts[idx]
                #if the point does not belong to the base shape, set its height (user-defined parameter)
                if pt not in lst_vertex: 
                    #v=context.bm.verts.new((pt.x, pt.y, height)) 
                    v=context.bm.verts.new((pt.x, pt.y, pt.z)) 
                    top_pts_list.add((pt.z, pt.y, pt.x))
                else: 
                    v=context.bm.verts.new((pt.x, pt.y, pt.z)) 
                v_lst.append(v)
            v_lst=tuple(v_lst)

            #add the face to the mesh 
            shape = createShape2d(v_lst)
            self.manager.resolve(shape)
        
        top_pts_list = list(top_pts_list)
        top_pts_list.sort(reverse = True)
        offset = top_height_z - top_pts_list[0][0]  
        for f in context.bm.faces: 
            for pt in f.verts: 
                if (pt.co.z, pt.co.y, pt.co.x) in top_pts_list: 
                    #pt.co.z = pt.co.z + offset
                    pass



class Corner:
    def __init__(self, vert, edge1, edge2, axis):
        self.vert = vert.co
        self._vert = vert
        self.edge1 = edge1
        self.edge2 = edge2
        # cross product between edge1 and edge1
        cross = edge1.vec.cross(edge2.vec)
        # To check if have a concave (>180) or convex angle (<180) between edge1 and edge2
        # we calculate dot product between cross and axis
        # If the dot product is positive, we have a convex angle (<180), otherwise concave (>180)
        dot = cross.dot(axis)
        self.convex = True if dot>0 else False
        # sine of the angle between -self.edge1.vec and self.edge2.vec
        sin = cross.length
        self.isLine = True if sin<zero and self.convex else False
        if not self.isLine:
            self.sin = sin if self.convex else -sin
            # cosine of the angle between -self.edge1.vec and self.edge2.vec
            self.cos = -(edge1.vec.dot(self.edge2.vec))
    
    def inset(self, d1, d2, translate=None, negate=False):
        vert = None
        if d1==0 and d2==0 and translate:
            vert = self.vert + translate
        else:
            if negate:
                d1 = -d1
                d2 = -d2
            vert = self.vert
            edge1 = self.edge1
            # extruded counterpart of self.vert
            vert = vert - d1*edge1.normal - (d2+d1*self.cos)/self.sin*edge1.vec
            if translate:
                vert = vert + translate
        if vert:
            self.vert = vert
            self._vert = context.bm.verts.new(vert)
    
    def updateForEvent(self, event, sequence):
        dt = event.t - self.t
        edge1 = self.edge1
        edge2 = self.edge2
        vert = self.vert - dt*edge1.distance*edge1.normal + dt*(-edge2.distance-edge1.distance*self.cos)/self.sin*edge1.vec
        if sequence.getVert:
            vert = sequence.getVert(vert, dt)
        self.vert = vert
        self.t = event.t
        self.event = event


class Edge:
    def __init__(self, vec, axis):
        self.vec = vec
        normal = vec.cross(axis)
        normal.normalize()
        self.normal = normal
    
    def length(self):
        return (self.corner2.vert - self.corner1.vert).length
    
    def ssInit(self):
        self.dirty = False
        self.edgeEvent = None
        self.leftVerts = [self.corner1._vert]
        self.rightVerts = [self.corner2._vert]


_id = 0
class Sequence:
    """
    Represents shrinking polygon for straight skeleton
    """
    def __init__(self, edge, numCorners, axis, getVert):
        global _id
        _id += 1
        self.axis = axis
        self.getVert = getVert
        self.id = _id
        # starting edge
        self.edge = edge
        self.numCorners = numCorners
        self.events = []
        self.eventIndex = 0
        
        # For edge with related pair of corners calculate intersection time.
        # If a straight skeleton is calculated for a roof,
        # the intersection time corresponds to a height
        while True:
            # add auxiliary lists to store vertices of the straight skeleton
            edge.ssInit()
            # also augment all corners
            edge.corner1.event = None
            edge.corner1.t = 0
            # calculate edge event
            self.addEdgeEvent(edge)
            edge = edge.corner2.edge2
            if edge == self.edge:
                break
       
    
    def process(self, changes):
        while True:
            # the current event
            event = self.events[self.eventIndex]
            event.resolve(self)
            self.eventIndex += 1
            if self.numCorners>2:
                edge = self.edge
                while True:
                    if edge.dirty:
                        if edge.edgeEvent:
                            self.removeEvent(edge.edgeEvent, self.eventIndex)
                        # update vert of edge.corner1 and edge.corner2
                        corner1 = edge.corner1
                        corner2 = edge.corner2
                        if corner1.event != event:
                            corner1.updateForEvent(event, self)
                        if corner2.event != event:
                            corner2.updateForEvent(event, self)
                        self.addEdgeEvent(edge, event.t)
                        edge.dirty = False
                    edge = edge.corner2.edge2
                    if edge == self.edge:
                        break
            else:
                break
    
    def addEdgeEvent(self, edge, offset=0):
            # left corner
            corner1 = edge.corner1
            # edge to the left from edge
            edge1 = corner1.edge1
            # right corner
            corner2 = edge.corner2
            # edge to the right from edge
            edge2 = corner2.edge2
            denominator = (edge1.distance+edge.distance*corner1.cos)/corner1.sin + (edge2.distance+edge.distance*corner2.cos)/corner2.sin
            # If denominator is equal to zero. it means that edge1 and edge2 are parallel.
            # If denominator is negative, it means that rays defined by edge1 and edge2 don't intersect
            if denominator>=zero:  
                t = edge.length()/denominator
                vert = corner1.vert - t*edge1.distance*edge1.normal + t*(-edge.distance-edge1.distance*corner1.cos)/corner1.sin*edge1.vec
                if self.getVert:
                    vert = self.getVert(vert, t)
                event = EventEdge(t+offset, edge, vert)
                edge.edgeEvent = event
                self.addEvent(event, self.eventIndex)
        
    
    def addEvent(self, event, lo=0):
        _lo = lo
        t = event.t
        events = self.events
        # implementing bisect.insort_right(..)
        hi = len(events)
        while lo < hi:
            mid = (lo+hi)//2
            if t < events[mid].t:
                hi = mid
            else:
                lo = mid+1
        # check if we have basically the same timestamp as events[lo]
        if lo < len(events):
            if events[lo].t-t <= timeTolerance:
                eventContainer = events[lo].append(event)
                if eventContainer:
                    events[lo] = eventContainer
                lo = -1
        # check if we have basically the same timestamp as events[lo-1]
        if lo>_lo:
            if t-events[lo-1].t <= timeTolerance:
                eventContainer = events[lo-1].append(event)
                if eventContainer:
                    events[lo-1] = eventContainer
                lo = -1
        if lo >= 0:
            events.insert(lo, event)
    
    def removeEvent(self, event, lo=0):
        container = event.container
        if container:
            container.remove(event)
            if container.numClusters == 0:
                self.removeEvent(container, lo)
        else:
            events = self.events
            hi = len(events)
            while True:
                mid = (lo+hi)//2
                if event == events[mid]:
                    break
                elif event.t < events[mid].t:
                    hi = mid
                else:
                    lo = mid+1
            del events[mid]


class EventContainer:
    def __init__(self, event1, event2):
        self.container = None
        self.events = []
        self.t = event1.t
        # init clusters
        self.clusters = []
        self.numClusters = 0
        self.append(event1)
        self.append(event2)
    
    def resolve(self, sequence):
        i = 0
        while i<self.numClusters:
            vert = self.clusters[i][0]
            self.clusters[i][0] = context.bm.verts.new(vert)
            i += 1
        for event in self.events:
            event.resolve(sequence, self.t, event.cluster[0])
    
    def append(self, event):
        self.clusterAssign(event)
        self.events.append(event)
        event.container = self
    
    def remove(self, event):
        self.events.remove(event)
        cluster = event.cluster
        cluster[1] -= 1
        if cluster[1]==0:
            # remove cluster
            self.clusters.remove(cluster)
            self.numClusters -= 1
    
    def clusterCreate(self, event):
        cluster = [event.vert, 1]
        event.cluster = cluster
        self.clusters.append(cluster)
        self.numClusters += 1
    
    def clusterAssign(self, event):
        clustered = False
        i = 0
        while i<self.numClusters:
            if self.clusterTest(event, i):
                self.clusterAdd(event, i)
                clustered = True
                break
            i += 1
        if not clustered:
            self.clusterCreate(event)
    
    def clusterTest(self, event, i):
        return (self.clusters[i][0]-event.vert).length_squared <= distanceTolerance2
    
    def clusterAdd(self, event, i):
        cluster = self.clusters[i]
        event.cluster = cluster
        cluster[1] += 1
    
    def __str__(self):
        return str(self.t) + ": " + str(len(self.events)) + " events"


class EventEdge:
    def __init__(self, t, edge, vert):
        # timestamp
        self.t = t
        self.edge = edge
        self.vert = vert
        self.container = None
    
    def resolve(self, sequence, t=None, vert=None):
        edge = self.edge
        if not t: t = self.t
        if not vert:
            vert = context.bm.verts.new(self.vert)
        # new corner
        corner = Corner(
            vert,
            edge.corner1.edge1,
            edge.corner2.edge2,
            sequence.axis
        )
        # set creation time (time offset) for the newly created corner
        corner.t = t
        # Mark the newly create edge, so we know that the edge has the correct vertex (vert).
        # We don't need to update vert later
        corner.event = self.container if self.container else self
        if edge.leftVerts[-1] != vert and edge.rightVerts[-1] != vert:
            edge.leftVerts.append(vert)
        sequence.numCorners -= 1
        # update the left neighbor
        edge1 = edge.corner1.edge1
        edge1.corner2 = corner
        edge1.dirty = True
        # update the right neighbor
        edge2 = edge.corner2.edge2
        edge2.corner1 = corner
        edge2.dirty = True
        
        if edge1.rightVerts[-1] != vert and edge1 != edge2:
            edge1.rightVerts.append(vert)
        if edge2.leftVerts[-1] != vert and edge1 != edge2:
            edge2.leftVerts.append(vert)
            
        # update the first edge of sequence
        if sequence.edge == edge:
            sequence.edge = edge2
    
    def append(self, event):
        return EventContainer(self, event)
    
    def __str__(self):
        return str(self.t) + " " + str(self.edge.corner1.vert)
