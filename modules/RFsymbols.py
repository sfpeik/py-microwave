from schemdraw.segments import *
import schemdraw.elements as elm
from numpy import sqrt, cos, sin, pi
       
class Transline(elm.Element2Term):
    def __init__(self, l, linecolor="black" ,*d, **kwargs):
        '''
        A transmission line as black bar 
        Behaves similar to a resistor twoport
        '''
        super().__init__(*d, **kwargs)
        self.segments.append(Segment( [(0, 0),(l, 0)],color="black",capstyle="square"))


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



class Antenna(elm.Element):
    ''' Antenna Fork Style '''
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        lead = 0.6
        h = 0.6
        w = 0.38
        self.segments.append(Segment([(0, 0), (0, lead), (-w, lead+h) ]))
        self.segments.append(Segment([        (0, lead), (w, lead+h)]))
        self.segments.append(Segment([        (0, lead), (0, lead+h)])) 
        self.elmparams['drop'] = (0, 0)
        self.elmparams['theta'] = 0
        self.anchors['start'] = (0, 0)
        self.anchors['center'] = (0, 0)
        self.anchors['end'] = (0, 0)

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
        llw = 6 * sqrt(size)
        self.segments.append(Segment([(0, 0), (size * 3, 0)], capstyle="butt", lw=6*size))
        self.segments.append(Segment([(0, -2*size), (size * 3, -2*size)], capstyle="butt", lw=llw))
        self.segments.append(Segment([(size * 0.7 , 0 * size), (size * 0.7, -2 * size)], lw=llw)) 
        self.segments.append(Segment([(size * 2.3 , 0 * size), (size * 2.3, -2 * size)], lw=llw)) 
        self.segments.append(Segment([(size * 0.8 , 0 * size), (size * 2.2, 0 * size)], lw=2*llw)) 
        self.segments.append(Segment([(size * 0.8 , -2 * size), (size * 2.2, -2 * size)], lw=2*llw)) 
        if showbox:
            self.segments.append(Segment([(0, 0.5 * size), (size * 3, 0.5 * size),(size * 3, -2.5 * size),(size * 0 , -2.5 * size),(0, 0.5 * size)], lw=2))

        self.anchors['p1'] = (0, 0)
        self.anchors['p2'] = (size*3, 0)
        self.anchors['p3'] = (size*3, -size*2)
        self.anchors['p4'] = (0, -size*2)
        self.params['drop'] = (size * 3, 0)
        self.params['pick'] = (size * 3, 0)
        
            
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
        
        

################################################################################################################ 
############# GENERAL  #########################################################################################
################################################################################################################
 
 
class Wigglyline(elm.Element):

    def __init__(self, ofst=(0,0), scale= 1, *d, **kwargs):
        '''
        Anchors:
        '''
        super().__init__(*d, **kwargs)
        size=0.3
        for x in range(0,6,2):
            self.segments.append(SegmentBezier([((0+x)*size,0),((0.5+x)*size, size),((1+x)*size,0)], capstyle="butt"))
            self.segments.append(SegmentBezier([((1+x)*size,0),((1.5+x)*size,-size),((2+x)*size,0)], capstyle="butt"))
        self.segments.append(Segment([((2+4)*size,0), ((2+4+1)*size,0)], capstyle="butt", arrow = '->'))
        
class Reflection(elm.Element):
    
    def __init__(self, ofst=(0,0), scale= 1, *d, **kwargs):
        '''
        Anchors:
        '''
        super().__init__(*d, **kwargs)
        dx,dy = ofst
        size = 0.6 * scale
        radius = 1.0
        lw = 3
        self.segments.append(Segment([(dx,dy), (0.5*size+dx,dy)], capstyle="butt"))
        self.segments.append(SegmentArc((0.5*size+dx,-radius*size/2+dy), radius*size, radius*size, theta1= -90, theta2=90, arrow = '<-'))
            

class ImpArrow(elm.Element):

    def __init__(self, len = 2, ofst=(-0.5,-1), *d, **kwargs):
        super().__init__(*d, **kwargs)
        x = ofst[0]
        y = ofst[1]
        self.segments.append(Segment([[0+x, 0+y], [-1+x, 0+y], [-1+x, -len+y]],arrow="<-") )

        
 
################################################################################################################ 
############# TWOPORTS #########################################################################################
################################################################################################################
 
 
class TwoportTline(elm.ElementTwoport):

    def __init__(self,lengthlabel="\\lambda/4", impedancelabel = "$Z_0$", **kwargs):
        linelabel="$\\leftarrow\\;"+lengthlabel+"\\; \\rightarrow$\n"+impedancelabel
        self.linelabel = linelabel
        super().__init__(input_element=elm.Gap, output_element=elm.Gap,boxpady=0.5, width=3.5, terminal=True, **kwargs)

    def setup(self):
        super().setup() 
        self.add(Transline(l=1,length=2.6).at(self.input_component.start).right().label(self.linelabel,loc="bot"))
        self.add(Transline(l=1,length=2.6).at(self.input_component.end).right())
        
        self.params['pick'] = ((-2,-2))
        self.drop(self.anchors['out_p'])
         
class TwoportPi(elm.ElementTwoport):
    """
    Pi-Network consisting of three Impedances 
    """

    def __init__(self, z1_element=elm.ResistorIEC, z2_element=elm.ResistorIEC, z3_element=elm.ResistorIEC, z1_label = "", z2_label = "", z3_label = "", **kwargs):
        self.z1_element = z1_element
        self.z2_element = z2_element
        self.z3_element = z3_element
        self.z1_label = z1_label
        self.z2_label = z2_label
        self.z3_label = z3_label
        
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.5, width=3.5, **kwargs)

    def setup(self):
        super().setup() 
        self.add(elm.Dot().at(self.input_component.start))
        self.add(self.z3_element(l=2).label(self.z3_label,loc="bot").at(self.input_component.start).right().dot())
        self.add(elm.Line(l=2.5).at(self.input_component.end).right())
        self.add(self.z1_element(l=1.5).label(self.z1_label,loc="bot",ofst=(.2, 0.1)).at(self.input_component.start).down().dot())
        self.add(self.z2_element(l=1.5).label(self.z2_label,ofst=(.2, -0.1)).at(self.output_component.start).down().dot())
        self.drop(self.anchors['out_p'])

class TwoportTee(elm.ElementTwoport):
    """
        Tee-Network consisting of three Impedances 
    """

    def __init__(self, z1_element=elm.ResistorIEC, z2_element=elm.ResistorIEC, z3_element=elm.ResistorIEC, z1_label = "", z2_label = "", z3_label = "", **kwargs):
        self.z1_element = z1_element
        self.z2_element = z2_element
        self.z3_element = z3_element
        self.z1_label = z1_label
        self.z2_label = z2_label
        self.z3_label = z3_label
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=3.5, **kwargs)

    def setup(self):
        super().setup() 
        self.add(self.z1_element(l=1.25).label(self.z1_label,loc="bot").at(self.input_component.start).right())
        self.add(elm.Line(l=.25).right().dot())
        self.add(self.z3_element(l=1.5).label(self.z3_label,loc="bot",ofst=(.2, 0.1)).down().dot())
        self.add(self.z2_element(l=1.25).label(self.z2_label,loc="bot").at(self.output_component.start).left())
        self.add(elm.Line(l=.25))
        self.add(elm.Line(l=3.5).at(self.input_component.end).right())
        self.drop(self.anchors['out_p'])


class TwoportShunt(elm.ElementTwoport):

    def __init__(self, shunt_element=elm.Resistor, shuntlabel="", **kwargs):
        self.shunt_element = shunt_element
        self.shlabel = shuntlabel
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=2.2, **kwargs)
          
    def setup(self):
        super().setup() 
        self.add(elm.Line(l=0.9).at(self.input_component.start).right().dot())
        self.add(self.shunt_element(l=1.5,label=self.shlabel).down().dot())
        self.add(elm.Line(l=.75).at(self.output_component.start).left())
        self.add(elm.Line(l=2.5).at(self.input_component.end).right())
        self.drop(self.anchors['out_p'])


class TwoportTransistor(elm.ElementTwoport):

    def __init__(self, trans_element=elm.Bjt, trans_label="", **kwargs):
        self.trans_element = trans_element
        self.trans_label = trans_label
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=2.2, **kwargs)
          
    def setup(self):
        super().setup() 
        self.add(elm.Line(l=0.8).at(self.input_component.start).theta(-75))
        self.add(elm.Line(l=0.2).right())
        self.add(self.trans_element(circle=True).right().label(self.trans_label))
        self.add(elm.Line(l=0.55).left().at(self.output_component.start))
        self.add(elm.Line(l=2.55).left().at(self.output_component.end))
        self.drop(self.anchors['out_p'])


class TwoportSeries(elm.ElementTwoport):

    def __init__(self, series_element=elm.Resistor, serieslabel="", **kwargs):
        self.series_element = series_element
        self.selabel = serieslabel
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=2.2, **kwargs)
        
    def setup(self):
        super().setup() 
        #self.add(elm.Line(l=1.75).at(self.input_component.start).right().dot())
        self.add(self.series_element(l=1.75).label(self.selabel,loc="bot").at(self.input_component.start).right())
        self.add(elm.Line(l=1.75).at(self.input_component.end).right())
        self.drop(self.anchors['out_p'])
        
class TwoportBlock(elm.ElementTwoport):

    def __init__(self, series_element=elm.Resistor, blocklabel="$[Z]$", **kwargs):
        self.series_element = series_element
        self.blocklabel = blocklabel
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=2.2, **kwargs)
    
    def setup(self):
        super().setup()
        self.segments.append(
                SegmentText(pos=Point((.85, .75)) + self.input_component.end, label=self.blocklabel,
                            align=('center', 'center'), fontsize=24, rotation_global=False))
        self.drop(self.anchors['out_p'])

class TwoportTransformer(elm.ElementTwoport):

    def __init__(self, series_element=elm.Resistor, translabel="$1:n$", **kwargs):
        self.series_element = series_element
        self.translabel = translabel
        super().__init__(input_element=elm.Gap, output_element=elm.Gap, boxpady=0.5,  boxpadx=0.0, width=2.2, **kwargs)
    
    def setup(self):
        super().setup()
        ll = 0.46
        self.add(elm.Line(l=ll).at(self.input_component.end).right())
        self.add(elm.Line(l=0.15).up())
        self.add(elm.Transformer(t1=3,t2=3,core=False).right().label(self.translabel,fontsize=11))
        self.add(elm.Line(l=ll).at(self.input_component.start).right())
        self.add(elm.Line(l=0.15).down())
        self.add(elm.Line(l=ll).at(self.output_component.start).left())
        self.add(elm.Line(l=0.15).down())
        self.add(elm.Line(l=ll).at(self.output_component.end).left())
        self.add(elm.Line(l=0.15).up())
        self.drop(self.anchors['out_p'])

class TwoportPort(elm.ElementCompound):
    """
    Creates two dots as port connections 
    """

    def __init__(self, portlabel="",labelcolor="k", **kwargs):
        self.portlabel = portlabel
        self.labelcolor=labelcolor
        super().__init__( **kwargs)

    def setup(self):
        self.add(elm.Dot(open=True))
        self.add(elm.Gap(l=1.5).down().label(self.portlabel,color=self.labelcolor))
        self.add(elm.Dot(open=True))

###################################################################################################################
###################################################################################################################
###################################################################################################################


if __name__ == "__main__":
    import schemdraw as schem
    import schemdraw.dsp as dsp
    
    choice = 5
    
    d = schem.Drawing()
    if choice == 1:
	    d += Port()
	    d += elm.Line(l=4)
	    d += elm.Resistor().down()
	    d += elm.Ground()
	    d += Reflection(color="r",dx= -1,dy=-1.5).label(r"$\Gamma_2$",fontsize=22)
	    d.draw()
    elif choice == 2:    
        d += (rat :=Ratrace().anchor('sum'))
        d += Port().label("In1",'right').fill('skyblue').at(rat.in1).down()
        d += Port(theta=30).label("In2",'left').fill('skyblue').at(rat.in2)
        d += Port(theta=-30).label(r"$\Sigma$",'left').fill('skyblue').at(rat.sum)
        d += Port().label(r"$\Delta$",'left').fill('skyblue').at(rat.delta).up()
        d.draw()
    elif choice == 3:     
        d += ( cir1 := Circulator('cw').flip().anchor('p1').fill('peachpuff') )
        d += Port(direction='right').up().at(cir1.p3).label("Port2").fill('skyblue')
        d += elm.Line(l=1).at(cir1.p2)
        d += elm.Coax(l=1,radius= 0.2).label(r"$\lambda/4$")
        d += (wilk := Splitter(size=1.3).label("Wilkinson\nDivider\n",fontsize=10).fill('khaki'))
        d += elm.Line(l=1)
        d += Isolator().fill('lightgray')
        d += elm.Line(l=1)
        d += (cou := Coupler(size=1.3).label("Forward\nCoupler\n10 dB\nno Phase Shift", fontsize=10).fill('khaki'))
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
        d.draw()
        d.save('example.png',dpi=600)
    elif choice ==5: 
        elm.style(elm.STYLE_IEC)
        d += TwoportPort()
        d.move(0.6,0)
        d += TwoportSeries(elm.Inductor,serieslabel="L")
        d += TwoportShunt(elm.Diode,serieslabel="L",boxls=":",box="r").label("Detector")
        d += TwoportPi(z1_element=elm.Capacitor, z2_element=elm.Capacitor, z3_element=elm.Inductor,
                       z1_label="$C_1$"        , z2_label="$C_3$"        , z3_label="$L_2$")
        d += elm.CurrentTransactor(terminals=True, boxpady=0.5,boxfill="lightgray"); d.move(2.6,0)
        d += TwoportTee(z1_element=elm.Capacitor, z2_element=elm.Capacitor, z3_element=elm.Inductor,
                       z1_label="$C_1$"        , z2_label="$C_3$"        , z3_label="$L_2$")
        d += TwoportTline()
        d += TwoportTransistor(boxfill="LemonChiffon")
        d += TwoportTransformer(translabel="1:4") 
        d.move(-0.2,0)
        d += TwoportPort()
        d += d.draw()     
       

        


