from pro import *

height= param(3.13)






@rule
def Begin():
    extrude(height,
        keepOriginal=True
    )
 