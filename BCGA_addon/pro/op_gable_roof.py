from .base import Operator, ComplexOperator, context

def gable_roof(*args, **kwargs):
    return context.factory["GableRoof"](*args, **kwargs)


class GableRoof(ComplexOperator):
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
        #if we have more than 2 numerical values, the first one is height and the others are soffits 
        if numValues>2:
            soffits = []
            height = args[0]
            for i in range(1,numValues):
                soffits.append(args[i])

        #if we have less than 2 numerical values or only one than the first one is the height and the others are soffits 
        else:
            height = args[0]
            if numValues==2:
                soffits = (args[1],)
            elif self.soffitSize:
                soffits = (self.soffitSize,)
        self.soffits = soffits
        self.top_height = height
        self.roof_height=None 
        
