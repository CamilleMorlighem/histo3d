import math
import bpy, bmesh
from pro import context

from .shape import createShape2d, createRectangle

from .util import zero
from .internal import straightSkeletonOfPolygon

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
        """
        #normal code: 
        # create faces and BCGA shapes for the straight skeleton
        for edge in self.edges:
            edge.leftVerts.reverse()
            face = [edge.leftVerts.pop()]
            face += edge.rightVerts
            face += edge.leftVerts
            shape = createShape2d(face)
            print("face")
            print(face)
            self.manager.resolve(shape, edge.pitch)
        """
        #dutch roof from here 
        faces_gabble_roof=[]
        for edge in self.edges:
            edge.leftVerts.reverse()
            face = [edge.leftVerts.pop()]
            face += edge.rightVerts
            face += edge.leftVerts   
            faces_gabble_roof.append(face)   
        
        for face in faces_gabble_roof: 
            if len(faces)==3: 
                starting_triang_face=face 
                break 
        
        for pt in starting_triang_face: 
            print(pt.co)

        





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


class Roof(Polygon):
    def __init__(self, verts, axis, manager=None):
        super().__init__(verts, axis, manager)

    def roof(self, base_shape, height):
        manager = self.manager
        face=base_shape.face
    
        """
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
        """
        """
        data = bpy.data.meshes.new(name='Straight Skeleton')
        obj = bpy.data.objects.new('Straight Skeleton', data)
        obj.location = bpy.context.scene.cursor.location
        bpy.context.scene.collection.objects.link(obj)
        obj.select_set(True)
        bpy.context.view_layer.objects.active = obj
        """
        lst_vertex=[]
        for v in face.verts: 

            lst_vertex.append(v.co)
   
        verts, edges, faces = straightSkeletonOfPolygon(polygon_vertices=lst_vertex, height=height)

        for face_range in faces: 
            v_lst=[]
            for idx in face_range: 
                pt=verts[idx]
                if pt not in lst_vertex: 
                    v=context.bm.verts.new((pt.x, pt.y, height)) 
                else: 
                    v=context.bm.verts.new((pt.x, pt.y, pt.z)) 
                print(pt.z)
                v_lst.append(v)
            v_lst=tuple(v_lst)
            print("vertices of one roof ace")
            print(v_lst)
            shape = createShape2d(v_lst)
            self.manager.resolve(shape)
            

                #create a list with the vertices for each face than call create2dshape 
                #than pass return shape in manager.resolve 

                #remove the steps where u create mesh straigt skeleton 
                #remove argument mesh data from straightskel fction 

        #obj.matrix_world = obj.matrix_world@plane_matrix

        #bpy.ops.mesh.straight_skeleton()
    
    


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
        for event in self.events:
            print(event)
    
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
            print("camiii")
            print(events)
            print("the event wer looking for")
            print(event)
            for x in range(0, len(events)): 
                print(events[x])
            
            while True:
                mid = (lo+hi)//2
                print(mid)
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
