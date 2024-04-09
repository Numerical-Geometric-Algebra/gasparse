-----------------
Multivector Types
-----------------

Depending on the algebra or on the data that is going to be handled we can chose from three different multivector types. The available multivector types are:

.. code-block:: c

    // Use dtype='sparse' to chose this type
    typedef struct SparseMultivector{
        int *bitmap;
        ga_float *value;
        Py_ssize_t size;
    }SparseMultivector;

    // Use dtype='blades' to chose this type
    typedef struct BladesMultivector{
        SparseMultivector *data;
        Py_ssize_t *grade;
        Py_ssize_t size;
    }BladesMultivector;

    // Use dtype='dense' to chose this type
    typedef struct DenseMultivector{
        ga_float *value;
        Py_ssize_t size; 
    }DenseMultivector;


``DenseMultivector`` are indexed by bitmaps. The ``BladesMultivector`` are indexed by grades where for each grade we have a ``SparseMultivector`` structure. 
The ``SparseMultivector`` represent multivectors by ``value`` and by its corresponding ``bitmap``, bitmaps indicate the basis multivector that corresponds to that value. 
The bit `i` indicates if there is a vector in the basis multivector. For example considering ``b'101101'`` the corresponding basis blade is :math:`e_{1346} = e_1\wedge e_3\wedge e_4\wedge e_6`. 
Using this representation for bitmaps makes it very easy to compute bitmaps resultant from the geometric product as ``b1 ^ b2`` where ``^`` is the bitwise xor operator. Yet signs are much more difficult 
to compute since it involves determining how many times we need to swap basis vectors to get to the canonical basis and determining the signs resulted from squaring vectors.

For the sparse multivectors we used an approach that allocates memory for the entire basis of the multivector and uses the indices as bitmaps. thus when we get a bitmap as the result of an operation we 
simply take that as an index e.g. ``dense.bitmap[dense0->bitmap[i]] = dense0.bitmap[i]``. Then we ignore all the unvisited indices by iterating over the entire multivector and checking if the bitmap 
is equal to ``-1``.

.. note:: 

    Iterating over the output multivector is expensive. For an `n` dimensional geometric algebra we have to loop over :math:`2^n` values and bitmaps.

For the ``'large'`` computation mode we compute bitmaps `online` by using the ``xor`` operator on the input bitmaps, then we use the same approach as with the ``'generic'`` computation mode.

.. admonition:: Question?

    How can I avoid iterating over all the multivector basis elements of a multivector?


To understand how we compute products for each of the different types of multivectors we show some functions for the different types


