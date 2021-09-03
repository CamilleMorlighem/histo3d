from pro import *

height= param(7)
top_height= param(12)





@rule
def Begin():
    extrude(height,
        top>>complex_hip_roof(top_height
        ),
        keepOriginal=True
    )
  