from pro import *

height= param(3.84)
top_height= param(7.23)
roof_pitch1= param(90)
roof_pitch2= param(30)
edge_midpt= (17.5, -1.0727, 0)


@rule
def Begin():
    extrude(height,
        top>>hip_roof2(top_height, roof_pitch1, roof_pitch2, edgeMidPoint=edge_midpt
        ),
        keepOriginal=True
    )
 