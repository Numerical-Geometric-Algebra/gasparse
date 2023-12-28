#!/usr/bin/env python

import jinja2

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader,trim_blocks=True)

# TEMPLATE_FILE_C = "multivector_array.c.src"
TEMPLATE_FILE_C = "largemultivector.c.src"

template_c = templateEnv.get_template(TEMPLATE_FILE_C)
text_op_funcs_c = template_c.render()


# to save the results

# with open("multivector_array.c", "w") as fh:
#     fh.write(text_op_funcs_c)
with open("largemultivector.c", "w") as fh:
    fh.write(text_op_funcs_c)
