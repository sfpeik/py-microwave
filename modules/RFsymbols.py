from schemdraw.segments import *
import schemdraw.elements as elm
from numpy import sqrt

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
        self.segments.append(SegmentArc((0,0), 2*radius*0.5, 2*radius*0.5, lw = 2.5, arrow=direction, theta1= 25, theta2=220))
        self.anchors['p1'] = (-off2,off2)
        self.anchors['p2'] = (off2,off2)
        self.anchors['p3'] = (0, -radius*leadextendfactor)
        self.params['drop'] = (off2,off2)
        self.params['pick'] = (off2,off2)
        self.params['lblloc'] = 'lft'

class Isolator(elm.Element):
    def __init__(self, *d, **kwargs):
        '''

        Anchors:
            * in
            * out
        '''
        super().__init__(*d, **kwargs)
        self.segments.append(Segment([(0, 0 ), (0.5,0)]))
        self.segments.append(Segment([(0.5, 1 ), (2.5, 1 )]))
        self.segments.append(Segment([(0.5, -1), (2.5, -1)]))
        self.segments.append(Segment([(0.5, 1), (0.5, -1)]))
        self.segments.append(Segment([(2.5, 1), (2.5, -1)]))
        self.segments.append(Segment([(2.5, 0), (3, 0)]))
        self.segments.append(SegmentArrow((0.9, 0), (2.3, 0), lw=2.5, headwidth=0.4, headlength=0.4))
        self.anchors['in'] = (0, 0)
        self.anchors['out'] = (3, 0)
        self.params['drop'] = (3, 0)

class Port(elm.Element):
    def __init__(self, direction='left', *d, **kwargs):
        '''
        direction: 'left' or 'right'

        Anchors:
            * p
        '''
        super().__init__(*d, **kwargs)
        size = 0.4
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

class Wilkinson(elm.Element):
    def __init__(self, showbox = True, *d, **kwargs):
        '''
        Anchors:
            * in
            *out1
            *out2
        '''
        super().__init__(*d, **kwargs)
        size = 1
        llw = 6
        self.segments.append(Segment([(0, 0), (size*1,0)],lw=5))
        self.segments.append(Segment([(size * 1, 0), (size * 2, size)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 1, 0), (size * 2, -size)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 2, size), (size * 3, size)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 2, -size), (size * 3, -size)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 2, size), (size * 2, 0.4*size)], lw=3, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 2, -size), (size * 2, -0.4* size)], lw=3, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 1.8,  0.4 * size), (size * 2.2, 0.4 *  size)], lw=2))
        self.segments.append(Segment([(size * 1.8, -0.4 * size), (size * 2.2, 0.4 * -size)], lw=2))
        self.segments.append(Segment([(size * 1.8,  0.4 * size), (size * 1.8, -0.4 *  size)], lw=2))
        self.segments.append(Segment([(size * 2.2, 0.4 * size), (size * 2.2, -0.4 * size)], lw=2))
        if showbox:
            self.segments.append(Segment([(-0.1, 1.5 * size), (size * 3+0.1, 1.5 * size)], lw=1))
            self.segments.append(Segment([(-0.1, -1.5 * size), (size * 3+0.1, -1.5 * size)], lw=1))
            self.segments.append(Segment([(-0.1, 1.5 * size), (-0.1 , -1.5 * size)], lw=1))
            self.segments.append(Segment([(3.1, 1.5 * size), (3.1 , -1.5 * size)], lw=1))
        self.anchors['in'] = (0,0)
        self.anchors['out1'] = (size * 3, size)
        self.anchors['out2'] = (size * 3, -size)
        self.params['drop'] = (size * 3, size)
        self.params['pick'] = (size * 3, size)


class Coupler(elm.Element):
    def __init__(self, showbox = True, *d, **kwargs):
        '''
        Anchors:
            * p1
            * p2
            * p3
            * p4
        '''
        super().__init__(*d, **kwargs)
        size = 0.5
        llw = 6
        self.segments.append(Segment([(0, 0), (size * 3, 0)], lw=6*size))
        self.segments.append(Segment([(0, -2*size), (size * 3, -2*size)], lw=6*size))
        self.segments.append(SegmentArrow((size * 0.3, -0.3 * size), (size * 2.6, -1.7 * size), lw=2*size)) 
        self.segments.append(SegmentArrow((size * 0.3, -1.7 * size), (size * 2.6, -0.3 * size), lw=2*size))
        self.segments.append(SegmentArrow((size * 2.6, -1.7 * size), (size * 0.3, -0.3 * size), lw=2*size)) # Reverse directions
        self.segments.append(SegmentArrow((size * 2.6, -0.3 * size), (size * 0.3, -1.7 * size), lw=2*size))
        if showbox:
            self.segments.append(Segment([(size * 0, 0.5 * size), (size * 3, 0.5 * size)], lw=1))
            self.segments.append(Segment([(size * 0, -2.5 * size),(size * 3, -2.5 * size)], lw=1))
            self.segments.append(Segment([(size * 0, 0.5 * size), (size * 0 , -2.5 * size)], lw=1*size))
            self.segments.append(Segment([(size * 3, 0.5 * size),  (size * 3 , -2.5 * size)], lw=1*size))

        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (size*3, 0)
        self.anchors['p3'] = (size*3, -size*2)
        self.anchors['p4'] = (0, -size*2)
        self.params['drop'] = (size * 3, 0)
        self.params['pick'] = (size * 3, 0)
        
 
if __name__ == "__main__":
    import schemdraw as schem
    d = schem.Drawing()
    d += Port().label("Port1")
    d += elm.Line(l=1)
    d += ( cir1 := Circulator('cw').flip().anchor('p1') )
    d += Port(direction='right').up().at(cir1.p3).label("Port2")
    d += elm.Line(l=1).at(cir1.p2)
    d += elm.Coax(l=1).label("$\lambda/4$")
    d += (wilk := Wilkinson().label("Wilkinson\nDivider\n",fontsize=10) )
    d += Isolator().reverse()
    d += (cou := Coupler().label("Forward\nCoupler\n10 dB\nno Phase Shift", fontsize=10) )
    d.push()
    d += elm.Line(l=3).at(wilk.out2)
    d.pop()
    d += elm.Coax(l=1).label("$\lambda/4$")
    d += Isolator().label("Isolator", fontsize=10)
    d += Port().left().label("Port3", loc='right')
    d += Isolator().at(cou.p3).reverse()
    d += Port().left().label("Port4", loc='right')
    d.save('example.png',dpi=600)
