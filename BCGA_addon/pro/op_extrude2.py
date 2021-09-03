from .base import Operator, ComplexOperator, context
import pro

section = "section"
middle = "middle"
cap1 = "cap1"
cap2 = "cap2"
cap = "cap"

def extrude2(*args, **kwargs):
    return context.factory["Extrude2"](*args, **kwargs)

class Extrude2(ComplexOperator):
    def __init__(self, *args, **kwargs):
        # is the first coordinate (the coordinate along the reference edge) relative or absolute?
        self.relativeCoord1 = True
        # is the second coordinate (the coordinate perpendicular to the reference edge) relative or absolute?
        self.relativeCoord2 = False
        
        self.inheritMaterialAll = True
        self.inheritMaterialSection = False
        self.inheritMaterialCap = False
        # Note: the original shape is kept irrespective of self.keepOriginal if
        # self.original is set
        self.keepOriginal = False
        self.symmetric = True
        self.axis = pro.x
        # shall we the normal for each extruded section?
        self.flipNormal = False
        # a rule for extruded section
        self.section = None
        # a rule for the last section
        self.last = None
        # a rule for the first cap
        self.cap1 = None
        # a rule for the second cap
        self.cap2 = None
        # a rule for both caps
        self.cap = None
        # a rule fo the original shape
        self.original = None
        # a rule for the middle face on for the symmetric case (self.symmetric==True)
        self.middle = None
        
        # apply kwargs
        for k in kwargs:
            setattr(self, k, kwargs[k])
        # process *args
        numOperators = 0
        # a list for parts definitions
        parts = []
        # are we in the declaration of parts or in the declaration of rules?
        inParts = True
        numArgs = len(args)
        argIndex = 0
        while argIndex<numArgs:
            arg = args[argIndex]
            if inParts:
                # Append to parts a tuple with 3 elements (value, height, operator) if the part is define like
                # value, height>>operator(...), ...
                # Append to parts a tuple with 2 elements (value, height) if an operator for the part is not given:
                # value, height, ...
                if isinstance(arg, Operator):
                    inParts = False
                else:
                    nextArg = args[argIndex+1]
                    if isinstance(nextArg, Operator):
                        part = (arg, nextArg.value, nextArg)
                        if nextArg.count:
                            nextArg.count = False
                            numOperators += 1
                    else:
                        part = (arg, nextArg)
                    parts.append(part)
                    argIndex += 2
            else:
                setattr(self, arg.value, arg)
                if arg.count:
                    arg.count = False
                    numOperators += 1
                argIndex += 1
        # Update parts for self.symmetric = True and self.relativeCoord1 = True
        # The case self.symmetric = True and self.relativeCoord1 = False will be treated in self.execute(),
        # since we don't know right now the absolute shape size
        if self.symmetric and self.relativeCoord1:
            argIndex = len(parts)-1
            # Treat the special case when there is no segment crossing 0.5 and parallel to the reference axis
            if parts[-1][0]!=0.5:
                # add the middle part
                part = parts[-1]
                parts.append((1-part[0], part[1]))
            while argIndex>0:
                argIndex -= 1
                part = parts[argIndex]
                prevPart = parts[argIndex+1]
                part = (1-part[0], part[1]) if len(prevPart)==2 else (1-part[0], part[1], prevPart[2])
                parts.append(part)
        self.parts = parts
        super().__init__(numOperators)


def extrude2_arc(radius, numSegments):
    pass