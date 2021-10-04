from schemdraw.segments import *
import schemdraw.elements as elm

class Circulator(elm.Element):
    def __init__(self, direction='cw', *d, **kwargs):
        '''
        direction: 'cw' or 'ccw'

        Anchors:
            * p1
            * p2
            * p3
        '''
        super().__init__(*d, **kwargs)
        radius = 1.0
        leadextendfactor = 1.5
        off = radius/sqrt(2)
        off2 = radius*leadextendfactor/sqrt(2)
        self.segments.append(SegmentCircle((0, 0), radius))
        self.segments.append(Segment([(-off, off), (-off2, off2)]))
        self.segments.append(Segment([(off, off), (off2, off2)]))
        self.segments.append(Segment([(0, -radius), (0, -radius*leadextendfactor)]))
        self.segments.append(SegmentArc((0,0), 2*radius*0.5, 2*radius*0.5, arrow=direction, theta1= 25, theta2=220))
        self.anchors['p1'] = (-off2,off2)
        self.anchors['p2'] = (off2,off2)
        self.anchors['p3'] = (0, -radius*leadextendfactor)
        self.params['drop'] = [off2,off2]
        self.params['pick'] = [off2,off2]
        self.params['lblloc'] = 'lft'

class Port(elm.Element):
    def __init__(self, direction='left', *d, **kwargs):
        '''
        direction: 'left' or 'right'

        Anchors:
            * p
        '''
        super().__init__(*d, **kwargs)
        size = 0.6
        dir = 1
        if direction =='right': dir = -1
        self.segments.append(Segment([(0, 0), (-size*dir, -size*0.7)]))
        self.segments.append(Segment([(0, 0), (-size*dir, size*0.7)]))
        self.segments.append(Segment([(-size*dir, -size*0.7), (-size*dir, size*0.7)]))
        self.anchors['p'] = (0,0)
        if direction == 'left':
            self.params['lblloc'] = 'lft'
        if direction == 'right':
            self.params['lblloc'] = 'rgt'
            
            
