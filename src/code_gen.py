#!/usr/bin/env python

import jinja2

data_types = ['blades','sparse','dense']
types = ['']

templateLoader = jinja2.FileSystemLoader(searchpath="./")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE_C = "operator_functions.c.src"
TEMPLATE_FILE_H = "operator_functions.h.src"
template_c = templateEnv.get_template(TEMPLATE_FILE_C)
template_h = templateEnv.get_template(TEMPLATE_FILE_H)
text_op_funcs_c = template_c.render(data_type_list=data_types,type_list=types)
text_op_funcs_h = template_h.render(data_type_list=data_types,type_list=types,n_data_types=len(data_types))

print(text_op_funcs_h)

# to save the results
# with open("operator_functions.c", "w") as fh:
    # fh.write(output_from_parsed_template)
