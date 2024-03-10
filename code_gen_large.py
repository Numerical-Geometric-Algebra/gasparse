#!/usr/bin/env python

import jinja2

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader,trim_blocks=True)

TEMPLATE_FILE_C = "src/multivector_large.c.src"

template_c = templateEnv.get_template(TEMPLATE_FILE_C)
text_op_funcs_c = template_c.render()

# save the results
with open("src/multivector_large.c", "w") as fh:
    fh.write(text_op_funcs_c)
