from pro import *

height= param(3.84)
top_height= param(7.23)
edge_midpt= (10.5, 0, 0)




@rule
def Begin():
    extrude(height,
        top>>gable_roof(top_height, edgeMidPoint=edge_midpt
        ),
        keepOriginal=True
    )
 