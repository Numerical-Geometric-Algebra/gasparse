gasparse.mvarray
----------------


.. py:class:: gasparse.mvarray

   The class where multivector arrays are instanciated. This class does not provide object creation. To initialize a multivector array use the method
   ``gasparse.GA.mvarray`` in :doc:`gasparse_ga`.

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
        
        Returns the grade of the multivector, up to some floating point precision. Returns an integer between ``0`` and ``p+q+r`` if it is of unique grade, otherwise returns ``None``. 
   
   .. py:method:: dual(self)
        
        Return the dualized multivectors.
 
   .. py:method:: undual(self)
        
        Return the undualized multivectors.

   .. py:method:: sum(self)
        
        Adds all the elements of the multivector array. Returns a 0-dimensional multivector array (aka multivector).

   .. py:method:: prod(self)

        Takes the geometric product of all the elements of the multivector array. Returns a 0-dimensional multivector array.

   .. py:method:: outer_prod(self)

        Takes the outer product of all the elements of the multivector array. Returns a 0-dimensional multivector array.

   .. py:method:: inner_prod(self)

        Takes the inner product of all the elements of the multivector array. Returns a 0-dimensional multivector array.
        **Carefull with precedence:** The product is taken by small indices first and the lower indices on the left. 
        e.g the inner product of an array of three multivectors computes ``((self[1]|self[2])|self[3])``

   .. py:method:: regressive_prod(self)

        Takes the regressive product of all the elements of the multivector array. Returns a 0-dimensional multivector array.
        **Carefull with precedence:** The product is taken by small indices first and the lower indices on the left. 
        e.g the regressive product of an array of three multivectors computes ``((self[1]&self[2])&self[3])``
 
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


The multivector array object ``gasparse.mvarray`` is described by understanding the core operations that can be done and what is dispatched when invoking operators like ``*``, ``+`` and ``^``. 
The operators invoke element wise operations on the multivectors of each argument. There are some restrictions on the use of these element wise operators.
The multivector arrays must have compatible metrics and their shape must be equal. For the metric to be compatible the metric arrays must overlap. 
Multivector arrays of different types or of different geometric algebras can also be added together, where they are cast to the biggest algebra.
To make implementation easier we do not make any distinction between multivectors and multivector arrays, other then the shape of the multivector arrays.
In fact multivectors are 0-dimensional multivector arrays.

Some operations are available only for multivector arrays of the ``scalar`` type. Those operations are:

- ``mvarray.cos(self)``
- ``mvarray.sin(self)``
- ``mvarray.cosh(self)``
- ``mvarray.sinh(self)``
- ``mvarray.exp(self)``
- ``mvarray.sqrt(self)``
- ``mvarray.sign(self)``

Note that ``mvarray.sign`` returns ``1`` when the input is ``0``.
