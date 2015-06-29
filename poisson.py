# -*- coding: utf-8 -*-
"""
The MIT License (MIT)

Copyright (c) 2015 Leonardo Kewitz
Copyright (c) 2015 Marcelo Vanti

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

from lib import *

nodes, elements = parsemesh('./res/rele.msh')
nn = len(nodes)
ne = len(elements)

mi0, eps0 = (4*pi*1E-7), 8.854187817E-12

# Configura materiais
# Material padrão Ar em todos os elementos.
for e in [es for es in elements]:
    e.mat = 1./(1.0*mi0)
# Muda o material do relè.
materiais = {
    24: 1./(1000*mi0),
    26: 1./(1000*mi0)
}
for superficie in materiais:
    valor = materiais[superficie]
    for elemento in elementsOnSurface(elements, superficie):
        elemento.mat = valor
# Adiciona densidade de carga no enrolamento.
for e in elementsOnSurface(elements, 22):
    e.f = 2000000

# Variáveis
Q = zeros(nn)
V = zeros(nn)
SS = zeros([nn, nn])
gradN = matrix('-1 1 0; -1 0 1')
N = array([1-1/3.-1/3., 1/3., 1/3.])

# Integração dos Elementos 2D.
for e in [es for es in elements if es.dim is 2]:
    n1, n2, n3 = [nodes[ni] for ni in e.nodes]
    indexes = [n1.i, n2.i, n3.i]
    # Calcula Jacobiano.
    J = matrix(zeros([2, 2]))
    J[0, :] = n2.x - n1.x, n2.y - n1.y
    J[1, :] = n3.x - n1.x, n3.y - n1.y
    # Inversa do Jacobiano.
    IJ = J.getI()
    # Gadiente N * Jacobiano.
    JGrad = IJ*gradN
    # Determinante do Jacobiano.
    detJ = linalg.det(J)
    # Integração final do elemento.
    intE = JGrad.getT()*JGrad*detJ*.5*e.mat
    # Right-side da equação.
    NjF = N*detJ*.5*e.f
    # Preenche matriz Q.
    for s in enumerate(indexes):
        Q[s[1]] += NjF[s[0]]
    # Preenche matriz de contribuição SS.
    for s in [(i, j) for i in enumerate(indexes) for j in enumerate(indexes)]:
        n1, n2 = s
        SS[n1[1], n2[1]] += intE[n1[0], n2[0]]

# Seta condições de contorno
contornos = {
    1: 0.0,
    2: 0.0,
    3: 0.0,
    4: 0.0,
    5: 0.0,
    6: 0.0
}
for linha in contornos:
    valor = contornos[linha]
    nos = nodesOnLine(elements, linha)
    for no in nos:
        Q[no] = valor
        SS[no, :] = 0.0
        SS[no, no] = 1.0

# Solve
V = linalg.solve(SS, Q)

# Plot
tri = triangulate(nodes)
plt.tricontourf(tri, V, 15, cmap=plt.cm.rainbow)
plt.colorbar()
plt.show()
