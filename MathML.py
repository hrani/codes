#Here I am passing expression to t and equivalent mathML expression are written out

import libsbml

t = "x1+2*x2"

t1 = libsbml.parseL3Formula(t)

t2 = libsbml.writeMathMLToString(t1)

print(t2)
'''
new_mathml = libsbml.parseL3Formula(t1)
new_string = libsbml.formulaToString(new_mathml)

print(new_mathml)
'''
