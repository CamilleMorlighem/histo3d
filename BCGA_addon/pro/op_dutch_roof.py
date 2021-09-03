from .base import Operator, ComplexOperator, context

def dutch_roof(*args, **kwargs):
    return context.factory["DutchRoof"](*args, **kwargs)


class DutchRoof(ComplexOperator):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.soffitSize = None
        self.fasciaSize = None
        self.n_steps = None 
        self.step_height = None 
        self.facade_thickness = None 
        self.edge_midpt = None 
        if "numberOfSteps" in kwargs: 
            self.n_steps = kwargs["numberOfSteps"]
        if "stepHeight" in kwargs: 
            self.step_height = kwargs["stepHeight"]
        if "fasciaSize" in kwargs:
            self.fasciaSize = kwargs["fasciaSize"]
        if "facadeThickness" in kwargs:
            self.facade_thickness = kwargs["facadeThickness"]
        if "edgeMidPoint" in kwargs: 
            self.edge_midpt = kwargs["edgeMidPoint"]

        self.face = None
        self.soffit = None
        self.fascia = None
        self.pitches = None 
        # find all definitions of operator and how many numerical values we have
        numOperators = 0
        numValues = 0
        i = len(args) - 1
        while i>=0:
            arg = args[i]
            if isinstance(arg, Operator):
                if isinstance(arg.value, str):
                    # rule or operator
                    setattr(self, arg.value, arg)
                elif numValues==0:
                    # roof definition (pitch or inset)
                    numValues = i + 1
                if arg.count:
                    arg.count = False
                    numOperators += 1
            elif numValues==0:
                # roof definition (pitch or inset)
                numValues = i + 1
            i -= 1
        self.numValues = numValues
        super().__init__(numOperators)
    
    def init(self, numEdges):
        args = self.args
        numValues = self.numValues
        soffits = None
        #if we have more than 2 numerical values, the first one is height, second is pitch and the others are soffits 
        if numValues>2:
            soffits = []
            height = args[0]
            pitches = [args[1]]
            for i in range(2,numValues):
                soffits.append(args[i])

            #if not soffits and self.soffitSize:
                #soffits = (self.soffitSize,)

        #if we have 2 numerical values than they are the height and the pitch 
        else:
            height = args[0]
            pitches=[args[1]]
            if self.soffitSize:
                soffits = (self.soffitSize,)
        self.soffits = soffits
        self.top_height = height
        self.roof_height=None 
        self.pitches=pitches 
        
