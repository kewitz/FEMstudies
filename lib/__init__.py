# -*- coding: utf-8 -*-
"""
The MIT License (MIT)

Copyright (c) 2015 Leonardo Kewitz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from numpy import *
from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import matplotlib.tri as tri


class Node(object):
    def __init__(self, args):
        self.i = int(args[0])-1
        self.x, self.y = float(args[1]), float(args[2])


class Element(object):
    def __init__(self, args):
        self.i = int(args[0])
        i, typ, ntags = args[:3]
        self.i, self.dim = int(i)-1, int(typ)
        self.tags = [int(a) for a in args[3:3+int(ntags)]]
        self.nodes = [int(a)-1 for a in args[3+int(ntags):]]


def readMesh(filename):
    """
    Processa um arquivo '.msh' do GMesh e retorna o conjunto (nós, elementos).

    >>> nodes, elements = readMesh('mesh.msh')
    """
    with open(filename) as f:
        x = f.read()
    ns, ne = x.find('$Nodes\n'), x.find('$EndNodes\n')
    nodes = [Node(n.split()) for n in x[ns+7:ne].split('\n')[1:-1]]
    es, ee = x.find('$Elements\n'), x.find('$EndElements\n')
    elements = [Element(e.split()) for e in x[es+10:ee].split('\n')[1:-1]]
    return nodes, elements


def nodesOnLine(*args):
    """
    Retorna índice dos nós existente em um elemento de linha.

    >>> boundary = nodesOnLine(elements, 0, 1, 2, 3)
    """
    elements = args[0]
    lines = args[1:]
    indexes = [ns for line in lines for es in elements
               if (es.dim is 1 and line in es.tags)
               for ns in es.nodes]
    return list(set(indexes))


def elementsOnSurface(*args):
    """
    Retorna os elementos existente em uma superfície.

    >>> sources = elementsOnSurface(elements, 10, 20)
    """
    elements = args[0]
    surfaces = args[1:]
    els = [es for surface in surfaces for es in elements
           if (es.dim is 2 and surface in es.tags)]
    return list(set(els))


def triangulate(nodes):
    points = array([(n.x, n.y) for n in nodes])
    x = points[:, 0]
    y = points[:, 1]
    return tri.Triangulation(x, y)
