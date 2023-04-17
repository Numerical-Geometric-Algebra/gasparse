


typedef struct _gen0_DenseMultivector{
    ga_float value[8]
}gen0_DenseMultivector;


PyObject *gen0_dense_multivector_geometric_product(PyMultivectorObject *data0, PyMultivectorObject *data1){
    gen0_DenseMultivector *pdense0 = data0->data;
    gen0_DenseMultivector *pdense1 = data1->data;
    gen0_DenseMultivector *pdense = (gen0_DenseMultivector*)PyMem_RawMalloc(sizeof(gen0_DenseMultivector));
    if(!pdense0 || !pdense1 || !pdense){
        PyMem_RawFree(pdense);
        return NULL; // raise error
    }

    gen0_DenseMultivector dense0 = *pdense0;
    gen0_DenseMultivector dense1 = *pdense1;
    gen0_DenseMultivector dense;
    PyMultivectorObject *out = new_multivector(data0);

    dense.value[0] =
    +dense0.value[0]*dense1.value[0]
    +dense0.value[1]*dense1.value[1]
    +dense0.value[2]*dense1.value[2]
    -dense0.value[3]*dense1.value[3]
    +dense0.value[4]*dense1.value[4]
    -dense0.value[5]*dense1.value[5]
    -dense0.value[6]*dense1.value[6]
    -dense0.value[7]*dense1.value[7]
;
    dense.value[1] =
    +dense0.value[0]*dense1.value[1]
    +dense0.value[1]*dense1.value[0]
    -dense0.value[2]*dense1.value[3]
    +dense0.value[3]*dense1.value[2]
    -dense0.value[4]*dense1.value[5]
    +dense0.value[5]*dense1.value[4]
    -dense0.value[6]*dense1.value[7]
    -dense0.value[7]*dense1.value[6]
;
    dense.value[2] =
    +dense0.value[0]*dense1.value[2]
    +dense0.value[1]*dense1.value[3]
    +dense0.value[2]*dense1.value[0]
    -dense0.value[3]*dense1.value[1]
    -dense0.value[4]*dense1.value[6]
    +dense0.value[5]*dense1.value[7]
    +dense0.value[6]*dense1.value[4]
    +dense0.value[7]*dense1.value[5]
;
    dense.value[3] =
    +dense0.value[0]*dense1.value[3]
    +dense0.value[1]*dense1.value[2]
    -dense0.value[2]*dense1.value[1]
    +dense0.value[3]*dense1.value[0]
    +dense0.value[4]*dense1.value[7]
    -dense0.value[5]*dense1.value[6]
    +dense0.value[6]*dense1.value[5]
    +dense0.value[7]*dense1.value[4]
;
    dense.value[4] =
    +dense0.value[0]*dense1.value[4]
    +dense0.value[1]*dense1.value[5]
    +dense0.value[2]*dense1.value[6]
    -dense0.value[3]*dense1.value[7]
    +dense0.value[4]*dense1.value[0]
    -dense0.value[5]*dense1.value[1]
    -dense0.value[6]*dense1.value[2]
    -dense0.value[7]*dense1.value[3]
;
    dense.value[5] =
    +dense0.value[0]*dense1.value[5]
    +dense0.value[1]*dense1.value[4]
    -dense0.value[2]*dense1.value[7]
    +dense0.value[3]*dense1.value[6]
    -dense0.value[4]*dense1.value[1]
    +dense0.value[5]*dense1.value[0]
    -dense0.value[6]*dense1.value[3]
    -dense0.value[7]*dense1.value[2]
;
    dense.value[6] =
    +dense0.value[0]*dense1.value[6]
    +dense0.value[1]*dense1.value[7]
    +dense0.value[2]*dense1.value[4]
    -dense0.value[3]*dense1.value[5]
    -dense0.value[4]*dense1.value[2]
    +dense0.value[5]*dense1.value[3]
    +dense0.value[6]*dense1.value[0]
    +dense0.value[7]*dense1.value[1]
;
    dense.value[7] =
    +dense0.value[0]*dense1.value[7]
    +dense0.value[1]*dense1.value[6]
    -dense0.value[2]*dense1.value[5]
    +dense0.value[3]*dense1.value[4]
    +dense0.value[4]*dense1.value[3]
    -dense0.value[5]*dense1.value[2]
    +dense0.value[6]*dense1.value[1]
    +dense0.value[7]*dense1.value[0]
;

    *pdense = dense;
    out->data = pdense;
    return (PyObject*)out;
}



