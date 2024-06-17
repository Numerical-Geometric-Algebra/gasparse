gasparse.GA
-----------

.. py:class:: gasparse.GA(p=3,q=0,r=0,metric=[1,1,1],print_type=0,print_type_mv=1,compute_mode="generic")
   
   Initializes a Geometric Algebra by either using the signature (p,q,r) where we have p positive, q negative and r null squaring basis vectors
   or either by defining the metric for the basis vectors, the value for which the basis vectors square to. 
   
   Returns the Geometric Algebra.

   :param p: number of positive squaring basis vectors.
   :param q: number of negative squaring basis vectors.
   :param r: number of null squaring basis vectors.
   :type p: integer
   :type q: integer
   :type r: integer
   :param metric: ``metric[i]`` is the metric of the i'th basis vector. The values can only be either (-1,0,1).
   :type metric: list[-1,0,1]
   :param print_type_mv: The print type of multivector arrays
   :param print_type: The print type of the algebra
   :param compute_mode: ``{'generic','large','generated'}`` computational mode.
        
        - ``'generic'``: tables are generated and stored for the different products
        - ``'large'``: only the table for the geometric product is generated the other products are computed online using the bitmaps
        - ``'generated'``: Algebras that were :code:`C` code generated for efficiency. Only few are available.
   :return: The Geometric Algebra.
   :rtype: gasparse.GA
   
   .. py:method:: size(self,grades=[],*args)

        Return the size of the algebra for the given grades, if no grades are give the size of the entire algebra is returned.
        Grades can be given as a variable number of arguments `x.size(1,2,3)` or as a list `x.size([1,2,3])`

        :param grades: A list of grades. Non-negative list of integers :code:`<= p+q+r`
        :type grades: list[int]
        :param args: A tuple of grades. Non-negative list of integers :code:`<= p+q+r`
        :type args: tuple[int] 
        :rtype: integer
        :return: The size of the algebra for all the grades in the grade list.
 

   .. py:method:: basis(self,grades=[],*args)

        Returns a dictionary for the basis blades of the algebra for the given grades, if no grades are give all the basis blades are returned.
        Grades can be given as a variable number of arguments or as a list.
        Use this to populate the python with the basis blades variables. e.g. ``locals().update(gasparse.GA(3).basis())``

        :param grades: A list of grades. Non-negative list of integers :code:`<= p+q+r`
        :type grades: list[int]
        :param args: A tuple of grades. Non-negative list of integers :code:`<= p+q+r`
        :type args: tuple[int]
        :rtype: dictonary{str:gasparse.mvarray}
        :return: Dictionary with the basis blades as values

   .. py:method:: set_precision(self,precision)

        Sets the precision for the removal of relatively small multivectors. After elementary operations the multivector basis elements are removed if the condition 
        ``values[i]*precision <= max(abs(values))`` holds true. Where ``values`` are the values of the multivector. This is done elementwise for each multivector in the multivector array.
        If the user wants to disable this option setting ``precision`` to zero does the trick.

        :param precision: The precision for which we decide to discard homogeneous multivector coordinates.
        :type precision: float
        :rtype: None

   .. py:method:: mvarray(values,grades=[],basis=[],dtype='sparse')
 
        Creates a multivector array from nested lists. The innermost lists sizes must agree with the size of the specified basis or grades. When grades and basis are not specified the whole algebra is considered. 
        Basis can be expressed either as list of strings, e.g. ``basis=['e1','e12']``, or as lists of multivectors, that is 0-dimensional multivector arrays.
        
        Returns a multivector array.
        
        :param values: The values of the basis elements.
        :type values: list[list] or list[float]
        :param basis: The basis blades associated with the values
        :type basis: list[gasparse.GA.mvarray] or list[str]
        :param dtype: ``{'sparse','dense','blades','scalar'}``, optional. A string indicating the type.
        
             * ``'sparse'`` basis elements are represented by arrays of values and corresponding bitmaps.
             * ``'dense'`` basis elements are represented by arrays of values the indices are its corresponding bitmaps.
             * ``'blades'`` basis elements are organized by grades where for each grade we have a ``'sparse'`` representation.
             * ``'scalar'`` the basis elements are one value which expresses the scalars of any geometric algebra. 
        
        
        :type dtype: string
        :return: The multivector array.
        :rtype: gasparse.mvarray


To initialize a geometric algebra you can use the ``gasparse.GA(p,q,r)`` function. Multivectors arrays are allways associated with some Geometric Algebra. 
The output of the ``gasparse.GA(p,q,r)`` function can then be used to create multivector arrays with the function ``gasparse.GA.mvarray``.
