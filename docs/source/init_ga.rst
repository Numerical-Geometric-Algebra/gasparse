Initializing Geometric Algebras and Creating Multivectors
---------------------------------------------------------

To initialize a geometric algebra you can use the ``gasparse.GA(p,q,r)`` function:

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


Multivectors arrays are allways associated with some Geometric Algebra. The output of the ``gasparse.GA(p,q,r)`` function can then be used to create multivector arrays with the function:


.. py:method:: gasparse.GA.mvarray(values,grades=[],basis=[],dtype='sparse')
   
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

The multivector array object ``gasparse.mvarray`` is described by understanding the core operations that can be done and what is dispatched when invoking operators like ``*``, ``+`` and ``^``. 
The operators invoke element wise operations on the multivectors of each argument. There are some restrictions on the use of these element wise operators. The multivector arrays must have compatible metrics and their shape must be equal. For the metric to be compatible the metric arrays must overlap. Multivector arrays of different types or of different geometric algebras can also be added together, where they are cast to the biggest algebra.
To make implementation easier we do not make any distinction between multivectors and multivector arrays, other then the shape of the multivector arrays. In fact multivectors are 0-dimensional multivector arrays.



.. py:class:: gasparse.mvarray

   The class where multivector arrays are instanciated. This class does not provide object creation. To initialize multivector use the method gasparse.GA.mvarray.

   .. py:method:: __add__(self,other)
      
      Adds the multivector array self with the multivector array other element wise. This function overloads the ``+`` operator.
     
      :type self: gasparse.mvarray
      :type other: gasparse.mvarray

      :return: A multivector array
      :rtype: gasparse.mvarray

   .. py:method:: __sub__(self,other)
        
        Subtracts other from self. Overloads the ``-`` operator.

   .. py:method:: __mul__(self,other)

        Computes the geometric product between two multivector arrays. Overloads the ``*`` operator. 

   .. py:method:: __xor__(self,other)

        Computes the outer product between self and other. Overloads the ``^`` operator.

   .. py:method:: __or__(self,other)
         
        Computes the inner product between self and other. Overloads the ``|`` operator.

   .. py:method:: __and__(self,other)
         
        Computes the regressive product between self and other. Overloads the ``&`` operator.

   .. py:method:: __truediv__(self,other)
         
        Divides ``self`` by ``other``. ``other`` must be of scalar type, float or int.  Overloads the ``/`` operator.
        
        :param other: A scalar type multivector array.
        :type other: gasparse.mvarray or float

   .. py:method:: __invert__(self)
         
        Applies grade reversion to the multivector array. Overloads the ``~`` operator.

   .. py:method:: grade(self)
        
        Returns the grade of the multivector, up to some float point precision.
   
   .. py:method:: dual(self)
        
        Return the dualized multivectors.
 
   .. py:method:: undual(self)
        
        Return the undualized multivectors.

   .. py:method:: sum(self)
        
        Adds all the elements of the multivector array. Returns a 0-dimensional multivector array (aka multivector).

   .. py:method:: prod(self)

        Takes the geometric product of all the elements of the multivector array. Return a multivector.
 
   .. py:method:: GA(self)

        Returns the algebra of a multivector array. 

        :return: The algebra associated with the multivector array.
        :rtype: gasparse.GA

   .. py:method:: type(self)

        Returns a string indicating the type of multivector array and the associated Geometric Algebra. 

   .. py:method:: tolist(self,grades=[],*args)

        Returns the multivector array as a list of floats and basis multivectors as ``gasparse.mvarray`` objects. 
        Grades can be passed as an argument to only output the values for those grades.


        :param grades: A list of grades. Non-negative list of integers :code:`<= p+q+r`
        :type grades: list[int]
        :param args: A tuple of grades. Non-negative list of integers :code:`<= p+q+r`
        :type args: tuple[int] 

        :return: A tuple with a nested list of values and a list of basis blades.
        :rtype: tuple(list[float], list[gasparse.mvarray])

   .. py:method:: cast(self,dtype)
        
        Casts the multivector array to the type ``dtype``. Can also be used to copy multivector arrays.

        :param dtype: ``{'sparse','dense','blades'}`` A string indicating the type.
        :type dtype: string
        :rtype: gasparse.mvarray
        :return: A multivector array of type ``dtype``

        
   .. py:classmethod:: concat(cls,list)

        Concatenates a list of multivectors. The shapes must agree otherwise it raises an error.
        Returns a multivector array. 


        :param list: A list of multivector arrays.
        :type list: list[gasparse.mvarray]
        :rtype: gasparse.mvarray


