from .base import Operator, ComplexOperator, context, countOperator

def join(neighbor, *args, **kwargs):
    return context.factory["Join"](neighbor, *args, **kwargs)

class Join(ComplexOperator):
    def __init__(self, neighbor, *args, **kwargs):
        self.neighbor = neighbor
        self.operator = None
        self.args = None
        numArgs = len(args)
        if numArgs: # if numArgs==1:
            operator = args[0]
            if isinstance(operator, Operator) and not hasattr(operator, "value"):
                self.operator = operator
            if numArgs>1:
                self.material = args[1]
        if not self.operator and numArgs:
            self.args = args
        self.kwargs = kwargs
        numOperators = 0
        for arg in args:
            if countOperator(arg):
                numOperators += 1
        super().__init__(numOperators)
