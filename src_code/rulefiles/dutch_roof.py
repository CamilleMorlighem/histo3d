from pro import *

height= param(3.84)
top_height= param(7.23)
roof_pitch= param(90)
step_height= param(0.4)
facade_thickness= param(0.5)
edge_midpt= (-8.75, 14.198, 0)

@rule
def Begin():
    extrude(height,
        top>>dutch_roof(top_height, roof_pitch, stepHeight=step_height, facadeThickness=facade_thickness, edgeMidPoint=edge_midpt
        ),
        keepOriginal=True
    )
