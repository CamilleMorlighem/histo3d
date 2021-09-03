from pro import *

height= param(7)
roof_height= param(3)
inset1= param(0.9)



   
@rule
def Begin():
    extrude(height,
        top>>inset(inset1,
            height=roof_height,
            keepOriginal=False 
        ), 
        keepOriginal=True 
    )
