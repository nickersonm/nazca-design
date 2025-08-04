#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module to easily create pin array's (a line of pins) with symbols and annotations.

Author: Ronald Broeke
Copyright (c) 2019-2023
"""
import nazca as nd
from nazca.logging import logger

class PinLine():
    """Class to define a set of pins on a line.

    The pinline can be used to place pins for pads or optical IO.
    The pinline is put in the layout inline with the pin in its put.

    * The pinline may shift along that pin using the balance attribute.

    A balance = 0.0 left aligns the pin line on its put location.
    A balance = 0.5 centers the pin line on its put location.
    A balance = 1.0 right aligns the pin line on its put location.
    The balance attribute may go outside the [0, 1] domain too.

    * Setting the length, pitch and number of pins N:

    The length of the pinline is set by either the length directly or by setting
    the pitch and number of pins (N). Similarly the pitch is set by either the pitch
    directly or by the length and N.

    * Pin names:

    Some examples:

    The pin name (and annotation and documentation text string) are by default
    constructed from a prefix attribute and an ordinal number,
    e.g. p0, p1, p2, ...

    The start (N0) and stepsize (dN) can be adapted
    e.g. for N0=3 and dN = 4 we get p3, p7, p11, ...

    For larger in pin counts the zero-fill parameter zfill can format the numbers more uniformly
    e.g. for zfill=3 we get p001, p002, p003, ...

    Change the order with reverse=True
    e.g. for N=10, N0=1, zfill=2, and reverse=True we get p10, p09, p08, ...

    Use prefix for change the default prefix 'p'
    e.g. for prefix='IO-', zfill=2, and N0=8 we get IO-08, IO-09, IO-10, ...

    It is possible to set a custom name per pin by providing attribute text_arr,
    a list of strings e.g. text_arr=['a', 5, hello'] and we get a, 5, hello.
    """
    def __init__(
            self,
            balance=0.5,
            length=None,
            N=None,
            dN=1,
            pitch=None,
            N0=0,
            flip=False,
            prefix='p',
            reverse=False,
            angle=0,
            symbol_cell=None,
            symbol_shape='arrow_full',
            symbol_scale=1.0,
            symbol_layer=500,
            text_move=[0, 0, -90],
            text_height= .5,
            text_show=False,
            text_arr=None,
            text_layer=None,
            text_align='cc',
            annotation=True,
            annotation_move=[-0.35, 0, 0],
            annotation_layer=None,
            instantiate=False,
            pintype=None,
            zfill=3,
            cellname='pinline',
        ):
        """Initialize Pinline object.

        The pin symbol can be set
        See method cell() for Args.

        Returns:
            None
        """
        self.pitch, self.N, self.length = self._configure(pitch, N, length)
        self.dN = dN                  # stepsize in counter
        self.N0 = N0                  # start value of ordinal counter
        self.balance = balance        # position where anchor the PinLine in a put 0=left, 0.5=middle, 1=right
        self.flip = flip              # flip pin direction in place
        self.prefix = prefix          # prefix of pin/annotation/text string
        self.reverse = reverse        # place pins in reverse order
        self.angle = angle            # rotate pins in place by angle
        self.symbol_cell = symbol_cell    # pin symbol cell
        self.symbol_shape = symbol_shape  # pin symbol polygon with a shape that will go into a cell
        self.symbol_scale = symbol_scale  # pin symbol scale
        self.symbol_layer = symbol_layer  # pin symbol layer
        self.text_move = text_move    # text displacement w.r.t. the pin
        self.text_height = text_height
        self.text_show = text_show    # show text if True
        self.text_arr = text_arr      # array to assign a custom text per pin.
        self.text_layer = text_layer  # text layer
        self.text_align = text_align  # text alignment string
        self.annotation = annotation
        self.annotation_move = annotation_move
        self.annotation_layer = annotation_layer
        self.instantiate = instantiate
        self.pintype = pintype
        self.zfill = zfill            # fill zeros to ordinal numbers
        self.cellname = cellname
        self.pin = []  # to fill for pin_iter (and skip auxiliary pins)
        if symbol_cell is None:
            self.symbol_cell = nd.bb_util.make_pincell(layer=symbol_layer, shape=symbol_shape)
        else:
            self.symbol_cell = symbol_cell


    def _configure(self, pitch, N, length):
        """Configure the PinLine dimensions.

        Note as least one of the arguments has to be None to avoid a conflict.

        Args:
            pitch (float): pitch
            N (int): number of pins
            length (float): length

        Returns:
            float, int, float: pitch, N, length
        """
        # pinline length and number parameters given
        # 3 parameters given:
        if pitch is not None and length is not None and N is not None:
            raise Exception("Only set one of the following combinations:"\
                "(length and pitch), (length and N), (pitch and N) not all:"\
                "length={}, pitch={}, N={}".format(length, pitch, N))
        # 2 given:
        elif length is not None and pitch is not None:
            N = int(length/pitch + 1)
            length0 = length
            length = pitch*(N-1) # lenght always adjusted to <= length set
            logger.debug(f'Adjust length from {length0} to {length} in pinline based on N and pitch')
        elif pitch is not None and N is not None:
            length = pitch * (N-1)
        elif N is not None and length is not None:
            if N > 1:
                pitch = length / (N-1)
            else:
                pitch = 0
        # 1 given:
        elif length is not None:
            N = self.N
            pitch = length / (N-1)
        elif pitch is not None:
            N = self.N
            length = pitch * (self.N-1)
        elif N is not None:
            pitch = self.pitch
            length = self.pitch * (N-1)
        # 0 given:
        else:
            length = self.length
            N = self.N
            pitch = self.pitch

        return pitch, N, length


    def cell(
            self,
            balance=None,
            N=None,
            N0=None,
            dN=None,
            length=None,
            pitch=None,
            prefix=None,
            zfill=None,
            reverse=None,
            flip=False,
            symbol_cell=None,
            symbol_scale=None,
            angle=None,
            angles=None,
            text_move=None,
            text_size=None,
            text_height=None,
            text_show=None,
            text_layer=None,
            text_align=None,
            text_arr=None,
            annotation=None,
            annotations=None,
            annotation_move=None,
            annotation_layer=None,
            pintype=None,
            instantiate=None,
            cellname=None,
    ):
        """Place pins.

        Create a line with equidistantly positioned pins. The default line
        lies on the positive x-axis, starting in x=0 and all pins looking in
        the positive y-direction. Note that either the length *or* pitch keyword
        can be used to define the pin distribution.

        The cell() parameters provided overule pinline settings for the specific
        cell placed.

        Args:
            balance (float): alignment of pin section.
                For example: 0.0 align on outer left of pinline,
                0.5 align on center, 1.0 align on outer right of pinline

            N (int): number of pins
            length (float): total length between outer pins in um
            pitch (float): pitch between pins in um

            prefix (str): optional annotation prefix before the pin number
            N0 (int): starting number
            reverse (bool): reverse numbering order
            flip (bool): flip the pin directions
            symbol_cell (Cell): pin symbol in layout

            angle (float): anlge of pins w.r.t. to default perpendicular.
            angles (list of float): custom angle per pin

            zfill (int): fill out number to zfill positions, e.g. 001, 002, etc for zfill=3
            annotation (bool): gds annotation output if True
            annotations (list of str): optional list to set the text of
                each pin explicitly. The list needs to have at least N elements
            annotation_move ((float, float, float)): offset of annotations w.r.t. the pin
            annotation_layer (str): layer for gds annotations

            text_show (bool): Show text output if True
            text_move ((float, float, float)): offset of text w.r.t. the pin
            text_layer (str): layer for text layer annotations
            text_align (str):
            instantiate (bool): instantiation value of cell with pins

        Returns:
            None
        """
        pitch, N, length = self._configure(pitch, N, length)

        if dN is None:
            dN = self.dN

        # other parameters
        if balance is None:
            balance = self.balance
#        if name is not None:
#            self.name = name
        if prefix is None:
            prefix = self.prefix
        if zfill is None:
            zfill = self.zfill

        if N0 is None:
            N0 = self.N0
        if reverse is None:
            reverse = self.reverse
        if flip is None:
            flip = self.flip
        if symbol_cell is None:
            symbol_cell = self.symbol_cell
        if symbol_scale is None:
            symbol_scale = self.symbol_scale
        if text_move is None:
            text_move = self.text_move
        if text_height is None:
            text_height = self.text_height
        if text_show is None:
            text_show = self.text_show
        if text_arr is None:
            text_arr = self.text_arr
        if text_layer is None:
            text_layer = self.text_layer
        if text_align is None:
            text_align = self.text_align
        if annotation is None:
            annotation = self.annotation
        if annotation_move is None:
            annotation_move = self.annotation_move
        if annotation_layer is None:
            annotation_layer = self.annotation_layer
        if instantiate is None:
            instantiate = self.instantiate
        if cellname is None:
            cellname = self.cellname
        if pintype is None:
            pintype = self.pintype
        if angle is None:
            angle = self.angle
        if flip:
            angle += -90
        else:
            angle += 90
        x0 = balance * length
        self.pos = []
        self.pin = []
        with nd.Cell(cellname, instantiate=instantiate, cnt=True) as C:
            for n in range(N):
                self.pos.append(length * n/(N-1))
                if not reverse:
                    i = n
                else:
                    i = N - n - 1
                if annotations is None:
                    anno = "{}{}".format(prefix, str(i*dN + N0).zfill(zfill))
                else:
                    anno = annotations[i]
                if text_arr is None:
                    text_str = anno
                else:
                    text_str = text_arr[n]
                p = nd.Pin(anno).put(-x0 + self.pos[-1], 0, angle)
                self.pin.append(p)
                symbol_cell.put(p, scale=symbol_scale)
                if text_show:
                    nd.text(
                        text_str,
                        height=text_height,
                        align=text_align,
                        layer=text_layer
                    ).put(p.move(*text_move))
                if annotation:
                    nd.Annotation(
                        text=anno,
                        layer=annotation_layer
                    ).put(p.move(*annotation_move))

            params = {
                'length': length,
                'pitch': pitch,
                'N': N,
                'N0': N0,
                'dN': dN,
                'reverse': reverse,
                'prefix': prefix,
                'zfill': zfill,
                'flip': flip,
            }
            C.properties['parameters'] = params
        return C


if __name__ == '__main__':
    # examples
    from pprint import pprint
    arrow = nd.bb_util.make_pincell(layer=2)
    arrow2 = nd.bb_util.make_pincell(layer=3, size=2.0)

    PL = PinLine()
    PL.cell(instantiate=False, balance=0).put(0)
    PL.cell(instantiate=False).put(0, 1)
    PL.cell(instantiate=False, balance=1).put(0, 2)
    PL.cell(instantiate=False, balance=0.25, reverse=True).put(0, 3)
    PL.cell(instantiate=False, balance=0.25, reverse=True, prefix='A').put(0, 4)
    PL.cell(instantiate=False, flip=True).put(0, 5)
    PL.cell(instantiate=False, flip=True, symbol_cell=arrow2).put(0, 6)
    PL.cell(instantiate=False, flip=True, pitch=4).put(0, 8)
    PL.cell(instantiate=False, pitch=4, N=4).put(0, 10)
    PL.cell(instantiate=False, pitch=4, N=4, dN=5).put(0, 11)
    PL.cell(instantiate=False, pitch=4, N=4, dN=5, N0=2).put(0, 12)
    PL.cell(instantiate=False, length=40).put(0, 13)

    PL = PinLine(symbol_shape='arrow_head', symbol_layer=6)
    PL.cell(instantiate=False, length=40, N=20).put(0, 14)
    PL.cell(instantiate=False, length=40, pitch=3).put(0, 15)
    PL.cell(instantiate=False, length=40, pitch=3, angle=25).put(0, 16)
    PL.cell(instantiate=False, length=40, pitch=3, angle=25, text_show=True, symbol_scale=4).put(0, 22)

    PL1 = PL.cell(
        symbol_cell=arrow,
        N=3,
        pitch=10,
        balance=0.0,
        prefix='H',
        annotations=['A', 'B', 'C'],
        text_show=True,
        text_height=1,
        text_move=(0, 0, 90),
        text_align='ct',
    )
    PL1.put(0, 30, 180)
    pprint(PL1.properties)


    # test pins in a cell template:
    L = 500
    N = 30
    buf = 25
    nd.Polygon(points=[(0, 0), (0, L), (L, L), (L, 0)], layer=3).put(0)


    # testmetal pin io
    pl1 = PL.cell(
        instantiate=False,
        length=L-4*buf,
        pitch=12,
        prefix='dcT',
        flip=True,
    ).put(0.5*L, L-buf, 0)

    N1 = pl1.properties['parameters']['N']
    pl2 = PL.cell(
        instantiate=False,
        length=L-4*buf,
        pitch=12,
        prefix='dcR',
        N0=N1,
        flip=True,
    ).put(L-buf, 0.5*L, -90)

    N2 = pl1.properties['parameters']['N']
    pl3 = PL.cell(
        instantiate=False,
        length=L-4*buf,
        pitch=12,
        N0=N1+N2,
        prefix='dcB',
        flip=True,
    ).put(0.5*L, buf, 180)


    # test waveguide edge-io with alternating angle
    pitch = 25
    pl4 = PL.cell(
        instantiate=False,
        pitch=pitch,
        N=6,
        dN=2,
        N0=0,
        prefix='ioW',
        zfill=2,
        flip=True,
        symbol_cell=arrow2,
        text_show=True,
        annotation=False,
        angle=0,
    ).put(0, 0.5*L, 90)
    pl5 = PL.cell(
        instantiate=False,
        pitch=pitch,
        N=6,
        dN=2,
        N0=1,
        prefix='ioW',
        zfill=2,
        flip=True,
        symbol_cell=arrow2,
        text_show=True,
        annotation=False,
        angle=7,
    ).put(0, 0.5*L+pitch*0.5, 90)

    nd.export_gds()
