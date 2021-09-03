from .base import Operator, ComplexOperator, context

def complex_hip_roof(*args, **kwargs):
    return context.factory["ComplexHipRoof"](*args, **kwargs)


class ComplexHipRoof(ComplexOperator):
    def __init__(self, *args, **kwargs):
        self.args = args
        self.soffitSize = None
        self.fasciaSize = None
        if "fasciaSize" in kwargs:
            self.fasciaSize = kwargs["fasciaSize"]
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
        #if we have more than 2 numerical values, the first one is height and the others are soffits 
        if numValues>2:
            soffits = []
            height = args[0]
            for i in range(1,numValues):
                soffits.append(args[i])

            #if not soffits and self.soffitSize:
                #soffits = (self.soffitSize,)

        #if we have less than 2 numerical values or only one than the first one is the height and the others are soffits 
        else:
            height = args[0]
            if numValues==2:
                soffits = (args[1],)
            elif self.soffitSize:
                soffits = (self.soffitSize,)
        self.soffits = soffits
        self.top_height = height
       