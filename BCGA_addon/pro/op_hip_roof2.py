from .base import Operator, ComplexOperator, context

def hip_roof2(*args, **kwargs):
    return context.factory["HipRoof2"](*args, **kwargs)


class HipRoof2(ComplexOperator):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.soffitSize = None
        self.fasciaSize = None
        self.edge_midpt = None 
        if "fasciaSize" in kwargs:
            self.fasciaSize = kwargs["fasciaSize"]
        if "edgeMidPoint" in kwargs: 
            self.edge_midpt = kwargs["edgeMidPoint"]
        self.face = None
        self.soffit = None
        self.fascia = None
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
        #when there are more than 2 numvalues, the first one is the height, the two nexts are the pitches and all the others are soffits 
        if numValues>3:
            soffits = []
            height = args[0]
            pitches = [args[1], args[2]]
            for i in range(3, numValues):
                soffits.append(args[i])
           
        #when there are only two values, the first one is the height and the second one is the pitch 
        else:
            height = args[0]
            if numValues==2: 
                pitches = [args[1], args[1]]
                if self.soffitSize:
                    soffits = (self.soffitSize,)
            if numValues == 3: 
                pitches = [args[1], args[2]]
                if self.soffitSize:
                    soffits = (self.soffitSize,)

        self.top_height = height 
        self.soffits = soffits
        self.pitches = pitches

