from schemdraw.segments import *
import schemdraw.elements as elm
from numpy import sqrt, cos, sin, pi


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
        radius = 0.5
        leadextendfactor = 1.5
        off = radius/sqrt(2)
        off2 = radius*leadextendfactor/sqrt(2)
        self.segments.append(SegmentCircle((0, 0), radius))
        self.segments.append(Segment([(-off, off), (-off2, off2)]))
        self.segments.append(Segment([(off, off), (off2, off2)]))
        self.segments.append(Segment([(0, -radius), (0, -radius*leadextendfactor)]))
        self.segments.append(SegmentArc((0,0), radius, radius, lw = 4*radius, theta1= -10, theta2=220))
        if direction == 'cw':
               self.segments.append(Segment([(radius*0.369,-0.15), (radius*0.369-0.07,-0.25)], arrow = '->'))
        if direction == 'ccw':
               self.segments.append(Segment([(-radius*0.369,-0.15), (-radius*0.369+0.07,-0.25)], arrow = '->'))
        self.anchors['p1'] = (-off2,off2)
        self.anchors['p2'] = (off2,off2)
        self.anchors['p3'] = (0, -radius*leadextendfactor)
        self.params['drop'] = (off2,off2)
        self.params['pick'] = (off2,off2)
        self.params['lblloc'] = 'lft'

class Isolator(elm.Element):
    def __init__(self, size= 1.0, *d, **kwargs):
        '''

        Anchors:
            * in
            * out
        '''
        super().__init__(*d, **kwargs)
        self.segments.append(Segment([(0, 0 ), (0,0.5), (1, 0.5 ), (1., -0.5 ), (0.0, -0.5), (0, 0)],lw=2))
        self.segments.append(Segment([(0.2, 0), (0.8, 0)], lw=2.5, arrow="->",arrowwidth=.2))
        self.anchors['in'] = (0, 0)
        self.anchors['out'] = (1, 0)
        self.params['drop'] = (1, 0)

class Port(elm.Element):
    '''
    RF Port Definition 
    
    '''
    _element_defaults = {
	'lead': True,
	'drop': (0, 0),
	'theta': 0,
	'lblloc': 'right'
    }
    def __init__(self, direction='left', *d, **kwargs):
        '''
        direction: 'left' or 'right'

        Anchors:
            * p
        '''
        super().__init__(*d, **kwargs)
        size = 0.4
        dir = 1
        if  (direction =='left'):
            dir = 1
        else:
            dir = -1 
        self.segments.append(Segment([(0, 0), (-size*dir, -size*0.7),(-size*dir, size*0.7),(0,0)]))
        self.anchors['p'] = (0,0)
        if dir==1:
            self.anchors['label'] = (-size * dir, 0)
        else:
            self.anchors['label'] = (0, 0)

class Splitter(elm.Element):
    def __init__(self, showbox = True, size=1 , couplabel="-3dB", *d, **kwargs):
        '''
        Anchors:
            * in
            *out1
            *out2
        '''
        super().__init__(*d, **kwargs)
        llw = 3*size
        self.segments.append(Segment([(0, 0), (size/2,0)],lw=llw, capstyle="butt",joinstyle= 'miter'))
        self.segments.append(Segment([(size * 0.5, 0),    (size * 1, size/2)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 0.5, 0),    (size * 1, -size/2)], lw=llw, joinstyle= 'miter'))
        self.segments.append(Segment([(size * 1, size/2), (size * 1.5, size/2)], lw=llw, capstyle="butt", joinstyle= 'miter'))
        self.segments.append(Segment([(size * 1, -size/2), (size * 1.5, -size/2)], lw=llw, capstyle="butt", joinstyle= 'miter'))
        ## Coupling Label
        self.segments.append(SegmentText((size*1,0),couplabel,fontsize = 9*sqrt(size)))
        if showbox:
            self.segments.append(Segment([(0, 0.75 * size), (size * 1.5, 0.75 * size),
                                          (size * 1.5, -0.75 * size),(0, -0.75 * size),
                                          (0, 0.75 * size)], lw=2))
        self.anchors['in'] = (0,0)
        self.anchors['out1'] = (size * 1.5, size/2)
        self.anchors['out2'] = (size * 1.5, -size/2)
        self.params['drop'] = (size * 1.5, size/2)
        self.params['pick'] = (size * 1.5, size/2)


class Ratrace(elm.Element):
    def __init__(self, showbox = True, size=1 , couplabel="Rat Race", *d, **kwargs):
        '''
        Anchors:
            * in
            *out1
            *out2
        '''
        super().__init__(*d, **kwargs)
        llw = 3*size
        radius = size
        an = 30*pi/180
        x,y = (-radius*cos(an), radius*sin(an))
        self.segments.append(SegmentCircle((0, 0), radius, lw=llw))
        self.segments.append(Segment([(0, radius), (0,radius + size/2)],lw=2*llw, capstyle="butt",joinstyle= 'miter'))
        self.segments.append(Segment([(0, -radius), (0,-radius - size/2)],lw=2*llw, capstyle="butt",joinstyle= 'miter'))
        self.segments.append(Segment([(x,  y), (x*(1+size/2), y*(1+size/2))], lw=2 * llw, capstyle="butt", joinstyle='miter'))
        self.segments.append(Segment([ (x,-y), (x*(1+size/2),-y*(1+size/2))], lw=2 * llw, capstyle="butt", joinstyle='miter'))
        ## Coupling Label
        self.segments.append(SegmentText((0,0),couplabel,fontsize = 9*sqrt(size)))

        self.anchors['sum'] = (x*(1+size/2), y*(1+size/2))
        self.anchors['delta'] = (0, -radius- size/2)
        self.anchors['in1'] = (0, radius+ size/2)
        self.anchors['in2'] = (x*(1+size/2),-y*(1+size/2))
        self.params['drop'] =  (-radius-size,radius)
        self.params['pick'] = (-radius-size,radius)
        
        
class Coupler(elm.Element):
    def __init__(self, showbox = True, size= 1, *d, **kwargs):
        '''
        Anchors:
            * p1
            * p2
            * p3
            * p4
        '''
        super().__init__(*d, **kwargs)
        size = 0.5 * size
        llw = 6
        self.segments.append(Segment([(0, 0), (size * 3, 0)], capstyle="butt", lw=6*size))
        self.segments.append(Segment([(0, -2*size), (size * 3, -2*size)], capstyle="butt", lw=6*size))
        self.segments.append(Segment([(size * 0.3, -0.3 * size), (size * 2.6, -1.7 * size)], arrow="<->", lw=2*size)) 
        self.segments.append(Segment([(size * 0.3, -1.7 * size), (size * 2.6, -0.3 * size)], arrow="<->", lw=2*size))
        if showbox:
            self.segments.append(Segment([(0, 0.5 * size), (size * 3, 0.5 * size),(size * 3, -2.5 * size),(size * 0 , -2.5 * size),(0, 0.5 * size)], lw=2))

        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (size*3, 0)
        self.anchors['p3'] = (size*3, -size*2)
        self.anchors['p4'] = (0, -size*2)
        self.params['drop'] = (size * 3, 0)
        self.params['pick'] = (size * 3, 0)
        
        
        
class Hybrid(elm.Element):
    def __init__(self, showbox = True, size= 1, *d, **kwargs):
        '''
        Anchors:
            * p1
            * p2
            * p3
            * p4
        '''
        super().__init__(*d, **kwargs)
        size = 0.5 * size
        llw = 6
        self.segments.append(Segment([(0, 0), (size * 3, 0)], capstyle="butt", lw=6*size))
        self.segments.append(Segment([(0, -2*size), (size * 3, -2*size)], capstyle="butt", lw=6*size))
        self.segments.append(Segment([(size * 0.7 , 0 * size), (size * 0.7, -2 * size)], lw=6*size)) 
        self.segments.append(Segment([(size * 2.3 , 0 * size), (size * 2.3, -2 * size)], lw=6*size)) 
        self.segments.append(Segment([(size * 0.8 , 0 * size), (size * 2.2, 0 * size)], lw=14*size)) 
        self.segments.append(Segment([(size * 0.8 , -2 * size), (size * 2.2, -2 * size)], lw=14*size)) 
        if showbox:
            self.segments.append(Segment([(0, 0.5 * size), (size * 3, 0.5 * size),(size * 3, -2.5 * size),(size * 0 , -2.5 * size),(0, 0.5 * size)], lw=2))

        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (size*3, 0)
        self.anchors['p3'] = (size*3, -size*2)
        self.anchors['p4'] = (0, -size*2)
        self.params['drop'] = (size * 3, 0)
        self.params['pick'] = (size * 3, 0)
        

class Reflection(elm.Element):
    
    def __init__(self, dx = 0, dy = 0 ,size= 1, *d, **kwargs):
        '''
        Anchors:
        '''
        super().__init__(*d, **kwargs)
        size = 0.8 * size
        radius = 1
        self.segments.append(Segment([(-0.5*size+dy, dx), (-0.5*size+dy, -0.5*size+dx)], capstyle="butt", lw=2))
        self.segments.append(SegmentArc((dy,dx), radius*size, radius*size, lw = 2, theta1= 0, theta2=180, arrow = '<-'))
        
        
class Attenuator(elm.Element):
    def __init__(self, variable = False, size= 1.0, *d, **kwargs):
        '''

        Anchors:
            * in
            * out
        '''
        super().__init__(*d, **kwargs)
        self.segments.append(Segment([(0, 0 ), (0,0.8), (1.6, 0.8 ), (1.6, -0.8 ), (0.0, -0.8), (0, 0)],lw=2))
        self.segments.append(Segment([(.8, 0.6), (0.6, 0.5), (1.0, 0.3), (.6, 0.1), (1.0, -0.1), (.6, -0.3), (1, -0.5), (.8, -0.6)], lw=1.5))
        if variable:
            self.segments.append(Segment([(0.2, -0.4), (1.4, 0.4)], lw=2, arrow="->",arrowwidth=.2))
        self.anchors['in'] = (0, 0)
        self.anchors['out'] = (1.6, 0)
        self.params['drop'] = (1.6, 0)
        
 
if __name__ == "__main__":
    import schemdraw as schem
    import schemdraw.dsp as dsp
    d = schem.Drawing()
    d += Port()
    d += elm.Line(l=4)
    d += elm.Resistor().down()
    d += elm.Ground()
    d += Reflection(color="r",dx= -1,dy=-1.5).label(r"$\Gamma_2$",fontsize=22)
    d.draw()
    exit(0)
    
    #d += (rat :=Ratrace().anchor('sum'))
    #d += Port().label("In1",'right').fill('skyblue').at(rat.in1).down()
    #d += Port(theta=30).label("In2").fill('skyblue').at(rat.in2)
    #d += Port(theta=-30).label(r"$\Sigma$").fill('skyblue').at(rat.sum)
    #d += Port().label(r"$\Delta$",'left').fill('skyblue').at(rat.delta).up()
    #d.draw()
    #plt.show()
    #exit(0)

    d += ( cir1 := Circulator('cw').flip().anchor('p1').fill('peachpuff') )
    d += Port(direction='right').up().at(cir1.p3).label("Port2").fill('skyblue')
    d += elm.Line(l=1).at(cir1.p2)
    d += elm.Coax(l=1,radius= 0.2).label(r"$\lambda/4$")
    d += (wilk := Splitter(size=1.3).label("Wilkinson\nDivider\n",fontsize=10).fill('Khaki'))
    d += elm.Line(l=1)
    d += Isolator().fill('lightgray')
    d += elm.Line(l=1)
    d += (cou := Coupler(size=1.3).label("Forward\nCoupler\n10 dB\nno Phase Shift", fontsize=10).fill('Khaki'))
    d.push()
    d += elm.Line(l=1).at(wilk.out2)
    d += dsp.Filter().fill('lavenderblush')
    d += elm.Line(l=1).right()
    d.pop()
    d += elm.Coax(l=1,radius=0.2).label(r"$\lambda/4$")
    d += dsp.Amp().label("LNA", fontsize=10).fill('linen')
    d += elm.Line(l=0.5)
    d += Port().left().label("Port3", loc='right').fill('skyblue')
    d += elm.Line(l=1).at(cou.p3)
    d += Isolator().reverse().fill('lightgray')
    d += elm.Line(l=2.5)
    d += Port().left().label("Port4", loc='right').fill('skyblue')
    d.save('example.png',dpi=600)


