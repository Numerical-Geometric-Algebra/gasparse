from geo_algebra import *

x = rdn_cga_vector_array(2)
y = rdn_cga_vector_array(2)
x_bar = x.sum()/100
B = x_bar ^ x
R = exp_trig(B)
R = exp_hyper(B)

mask = (x|x)(0) > (y|y)(0)