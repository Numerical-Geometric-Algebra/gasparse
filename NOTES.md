# sparse-multivectors

Creating a python library to do computations using sparse representation of multivectors.

One of the motivating things about writting a python `C API` is that I can dispatch python operations to `C` implemented code. Overriding python operators to dispatch `C` code is such a complicated and generalized approach to implement a GA that I do not know why have I ever started... Is it because I am an optimization freak that really likes programming in c or because I wanted to have my own approach to a GA software? 

Anyway I kinda wanted to copy numpy's approach, but instead of having just arrays of scalars could also have arrays of multivectors. The big difficulty in geometric algebra, in my opinion, is the computation of the sign wich results from the product of two basis blades. The computation of this quantity for big algebras can add a huge overhead. When signs are availaible through a table the cost of fetching that sign gives performance hits to the code. If there where a very performant way to compute the signs, it would make GA leap to a whole new level. The problem with the signs is that they take up space. The smallest scalar type we can use is a char which is `1` byte, and using `1` byte for a number which is either `-1`, `0` or `1` is rather memory ineffecient.

One of my implementation objectives was to consider sparse multivectors. A first approach was taken by considering bitmaps and value pairs for each multivector, bitmaps representing basis blades and values its corresponding scalar value. We initially compute binary tables for the signs and the bitmaps for the different products. The tables are then used to compute the different products. Since multivectors are usually grade sparse it would make sense to compute tables for each grade instead for the whole algebra. This way we do not need to load the whole table into memory(Maybe it does not anyway). It makes sense that somehow this should be more efficient.

Different multivector types. Due to efficiency considerations we decided to have multiple multivector types. Concretely for big algebras we consider sparse multivectors, For small algebras we can use dense multivectors or code generated types. By considering different multivector types we can also dispatch operations which are only available to certain types of multivectors. 


**Future Work:** Some multivectors have specific properties that makes them simple to do computations with. Thus we aim to provide special types of multivectors for which computations can be readily dispatched in the most efficient manner. Thus in future versions we aim to develop specialized `C` code for these types of operations, along with the different types of multivectors. 

Concretely multivectors $X$ which satisfy $XX^\dagger = \langle XX^\dagger\rangle$ are specially simple to use and do computations with specifically inversions can be computed as $X^{-1}=X/\langle XX^\dagger\rangle$. Thus if we specify the type of multivector which is guaranteed to have that property we can make negative powers of multivectors dispatch the operation $X^{-k} = {X}^k/\langle XX^\dagger\rangle^{k}$. Other types of multivector for which operations can be readily simplified are multivector that square to scalars $X^2=\langle X^2\rangle$ then logarithms and exponentials can be easily computed 
$\exp(X) = \cosh(X)+\sinh(X)$ where $\cosh(X) = \sum_{k=0}^\infty \frac{\sqrt{\langle X^2\rangle}^{2k}}{(2k)!} =\cosh(\sqrt{\langle X^2\rangle}) $ and $\sinh(X) = \hat{X} \sum_{k=0}^\infty \frac{\sqrt{\langle X^2\rangle}^{2k+1}}{(2k+1)!}=\hat{X}\sinh(\sqrt{\langle X^2\rangle})$ which depending on the sign of the square of the multivector gets dispatched to hyperbolic ($X^2>1$), trigometric ($X^2<1$) or parabolic ($X^2=0$). 

Computations of matrices of multilinear functions can be improved by writting specialized `C` code which given some multivector representation of the multilinear function can be optimized to compute $F_{ij} = f_i*F(f^j)$ where $f_i*f^j = \delta_{ij}$ are basis blades of the geometric algebra. 

Ternary product computations can also be optimized. In particular grade preserving ternary products can be optimized. In order to make this optimization we have to compute ternary tables for all combinations of grades of the algebra, in particular we construct a ternary table `map[grade0][grade1][grade2][grade_out]` for each grade. That multiplies `multivector0` with `multivector1` with `multivector2` of grades `grade0`, `grade1`, `grade2` respectively. The multiplication is done by numerical order from left to right `result <- multivector0*multivector1*multivector2`. For generated algebras (which are in general smaller) the efficiency improvements should be apparent, since the tables are smaller, and consequentely the number of function pointers associated with it are also small. In a generated algebra we would consider dispatching ternary operations to specific types of multivectors. The main problem with this approach for big algebras is that the ternary tables become very big which means that they will occupy more space on memory, also accessing this memory might become more expensive. To make better use of memory the graded ternary maps are only going to be computed for those specific grades, this way only when the graded table is to be used, is when it will be generated. *Question: What is the compuational cost difference of considering a graded ternary table VS a binary table, where in the later the product is done for all the grades of the multivector?* 

To provide efficient graded ternary product computations we will have to define two new types of multivectors which are going to be indexed by grade and by position. We will have to define graded dense multivectors and graded sparse multivectors. The ternary map will have to map (for some particular grades) will have to map three positions (the input multivectors) to a position (the output multivector). For example for dense multivectors we might consider the following: 
```c
typedef struct GradedDenseMultivector{// Use GradeMap to determine the size for each grade
    ga_float **values;
    Py_ssize_t *grades;
    Py_ssize_t size; // The number of different grades
}GradedDenseMultivector;
```

Signs for binary and ternary tables are expensive to compute and require memory storage accessing needed when computing products. In future work we want to construct an hash function which can determine signs in a very fast way. Improving on the binary and ternary maps approach by not needing to store big tables. Another improvement that can also be made is to compress signs to a single integer, thus a single integer would represent multiple signs. Using a trinary representation of the integers where `0` would represent `-1`, the value `1` would represent `0` and `2` would represent `1`. The n-th trit would represent the $n$-th sign.

### Repository TODO


## VERY IMPORTANT
**upload tar.gz to PyPI**

1. Create branch for multilinear algebra and for ternary operators.
1. Remove all content related with multilinear algebra from this branch

     
    1. Create multilinear repository (pushing the latest changes to it)
    2. Remove multilinear and ternary operations from the main branch.
	  2. Create copy of the readme file as NOTES.md
	  2. Write a new readme file
	  3. Clean code. Remove unnecessary code. Compile with optimizations enabled.
    3. Create a git tag for the main branch


## DONE

- Create a geometric algebra via p,q,r, or via a metric
```python
vga = gasparse.GA(3)
cga = gasparse.GA(4,1) 
cga = gasparse.GA(metric=[-1,1,1,1,1]) 
ga = gasparse.GA(p=2,q=3,r=4)
```
  - The order of the basis are:
    1. `p` positive basis vectors
    2. `q` negative basis vectors
    3. `r` null basis vectors

Note that the last two lines generate algebras with the same dimensions yet the symbols represent different basis, as such they are not compatible. Furthermore the last line is not compatible with the vga algebra, since the first basis element is negative. Note that the following will result in an **error**:
```python
vga = gasparse.GA(3)
cga = gasparse.GA(metric=[-1,1,1,1,1])
t = cga.multivector([1,2,3,4,5],grades=1)
y = vga.multivector([1,2,3],grades=1)
z = t + y
```

The user can chose between three computation modes `'generic'`, `'large'`, `'generated'`
To generate an algebra under these computation modes use
```python
ga = gasparsegen.GA(p,q,r,compute_mode=mode)
```
where `mode` is a string chosen from `'generic'`, `'large'`, `'generated'`.

  - The *generic* type stores the bitmaps in a bitmap array (these are generated when the algebra is initialized). 
    It uses the bitmaps array and the sign array to compute the different products
  - The *large* type generates only the signs at the initialization. It has a shorter initialization but computationaly slower (Good for large algebras)
  - The *generated* type are code generated geometric algebras (As of now the availables GA's are 3DVGA and 3DCGA)


- Initialize a basis
```python
basis = cga.basis()
locals().update(basis) # Update all of the basis blades variables into the environment
```
The basis will be available as the variables `e1`, `e2`,..., `e12`,...,`e12345`. 

- Initialize multivector via:
  - A given basis (multivector or string)
  ```python
  x = ga.multivector([1,23,4,4],[1.0,e1,e2,e12])
  y = ga.multivector([1,23,4,4],['e','e1','e2','e12'])
  ```
  
  - A given list of grades or a single grade
  ```python
  u = ga.multivector([-1,1,2,3],grades=[0,2])
  v = ga.multivector([1,2,3],grades=2)
  ```
  
  - Nothing (Considers all basis elements)
  ```python
  z = ga.multivector([1,2,3,4,5,6,7,8])
  ```
  - The user can also initialize the multivector in a some chosen type
  ```python
  x = ga.multivector([1,23,4,4],[1.0,e1,e2,e12],dtype='dense')
  ```
  - When the user wants all multivectors to be of the dense type then he can use the `ga.default` methods as: 
  ```python
  ga.default("dense")
  x = ga.multivector([1,23,4,4],[1.0,e1,e2,e12]) # This will be of type dense
  ```

Note that the values for the `ga.multivector` function can be either integer or float. (Everything gets cast to `ga_float`.)

- Function that returns sizes:
  - For each grade `ga.size(2)`
  - For multiple grades `ga.size([1,2])` or `ga.size(1,2)`
  - For the whole algebra (empty argument) `ga.size()`

- Zero grade projection now outputs objects of type float instead of `gasparse.multivector`
- The user can now use `x.grade()` to get the grade of a multivector if it has multiple grades then it returns `None`

- The function list accepts grades as input, can be a single integer or a list 
```python
values,basis = x.list([0,1]) # returns grades zero and one
values,basis = x.list(0,1) # returns grades zero and one
values,basis = x.list() # returns all elements
values,basis = x.list(1) # returns grade one with a list of the multivector basis
```
The first variable of the output are the values atributed to each basis, the second are the bitmaps that indicates the basis, the basis bitmaps can be chosen to be bitmaps or multivectors (in the later  case the basis is of the same dtype as `x`)

- Get the size of the algebra or the size of each grade
```python
ga.size() # gets the algebra size
ga.size(2) # Gets the size of the bivector basis
ga.size([1,2]) # Gets the sum of the size of  grade one and grade two
```

- All the operations defined in a geometric algebra are available, given two multivectors x and y in compatible algebras we can write
```python
x*y # geometric product
x|y # Inner product
x^y # wedge product
x&y # Regressive product
x+y # Addition
x-y # subtraction
a*x # scalar multiplication
~x # multivector reverse
x(g) # projection to grade g
x([g1,g2]) # projection to grades g1 and g2
x(g1,g2) # projection to grades g1 and g2
```

Grade projection to the zero grade returns a scalar type array.
```python
x = ga.multivector([[2,2,3,4],[5,6,7,8]],grades=[0,1])
print(x(0).Type()) # returns 'scalar'
```
Division by the scalar type is defined!

Note that the division `/` operator is only defined for division by scalars or scalar arrays since for general multivectors it is not simple to compute.
The following operations are valid
```python
 x/y(0)
 y(0)/12.34
 12.34/y(0)
 x/12.34
```
Note that `y(0)` will return a new type of multivector which is a scalar type multivector. Which can be handled as a scalar. Thus division by this quantity is possible. Then we can compute inverse of homogeneous multivectors arrays with 
```python
  x_inv = x/(x*~x)(0)
```

- Get the grade of a multivector
```python
x.grade()
```
Returns `None` if it is not of unique grade.

- Compute the dual and undual of a multivector
```python
y = x.dual()
z = y.undual()
```
It is the same as multiplying by the unit pseudoscalar when in a non-null metric. When we have a null basis element the dual is defined in a different way (expalain further...).

  - Atomic operations....
  - Code generated types...
  - Casting to another type/algebra...
    - The user can cast a multivector to some other type
    ```python
    x = ga.multivector([1,23,4,4],[1.0,e1,e2,e12],dtype='dense')
    y = x.cast("sparse")
    ```
    Note that casting from the dense type will not remove the basis elements which are zero. If the user wants to get rid of the zeros he can use `y *= 1` or just `y = 1*x.cast("sparse")`.

For the code generated algebras we have two types the `dense` and the `blades`. The key difference is in their internal representations. Both types are dense representations of a multivector. To use a specific type the user can call `y = x.cast('blades')` or using the equivalent for the other types such exemplified previously. The user cannot use the cast function to cast multivectors between different algebras in order to do that we will make available `ga.cast` which casts a multivector to the algebra ga.

 - Grade projection to the grade zero outputs a float, which can be used directly by numpy. That is, the output of 
 ```python
 x = ga.multivector(x_list)
 print(type(x(0)))
 ```
should be `float`.

Now the user can create arrays of multivectors for example it takes lists of lists
```python
ga = gasparse.GA(3)
basis = ga.basis()
locals().update(basis)
lst = [[[[1.1,2.2,3.3],[4.0,5.0,6.0]],
        [[7.0,8.0,9.0],[10.0,11.0,12.0]]],
       [[[1.0,2.0,3.0],[4.0,5.0,6.0]],
        [[7.0,8.0,9.0],[10.0,11.4,12.0]]],
       [[[1.0,2.0,3.0],[4.0,5.0,6.0]],
        [[7.0,8.0,9.0],[10.0,11.0,12.6]]],
       [[[1.0,2.0,3.0],[4.0,5.0,6.0]],
        [[7.0,8.0,9.0],[10.9999,11.0,12.0]]]]
        
y = ga.multivector(lst,basis=[e1,e2,e3])
```
Then to retrieve the exact same list the user can use
```python
lst,basis = y.list(1)
```
Note how the list function also outputs the basis.

 - To get information regarding a multivector the user can use
	```python
	x.GA() # Returns the geometric algebra of that element
	x.type() # Returns a string indicating the type

	```
	- `'sparse'`, `'dense'`, `'blades'` are the generic types.
	- `'dense0'`,`'dense1'`, ... ,`'blades0'`,`'blades1'`, ...  are the code generated types.
	- `'denselarge'`,`'sparselarge'`,`'bladeslarge'`  are the types that use the large computation mode.

	e.g.
	```python
	>>> import gasparse
	>>> ga = gasparse.GA(3)
	>>> x = ga.multivector([1,2,4],grades=1,dtype='sparse')
	>>> x.type()
	'GA(3).mvarray.sparse'
	```

  - Operations on a single multivector array
  ```python
  x = ga.multivector([[0,1],[2,1]],basis=['e1','e2'])
  z = x.prod() # Geometric Multiplies the element the first element is the leftmost element in the product
  z = x.outer_prod() # Outer Multiplies the element the first element is the leftmost element in the product
  z = x.sum() # Sums all the multivectors in the array
  ```
  Note that for the inner and regressive products we have precedence to consider, that is
  ```python
  z = a|b|c # This is the same as the next line
  z = (a|b)|c
  z = a|(b|c) # Which is not the same as this
  ```
  Since the inner and regressive product is more ambiguous the user has to be carefull when using it since the following relations hold
  ```python
  a|b|c|d <=> ((a|b)|c)|d
  a.inner_prod() <=> a[1]|a[2]|a[3]|a[4] <=> ((a[1]|a[2])|a[3])|a[4]
  ```
  If the user wants he can create scalar multivector arrays via the keyword `dtype='scalar'`, for example
  ```python
  x = ga.multivector([[1.234],[2.345]],grades=0,dtype='scalar')
  x = ga.multivector([[1.234],[2.345]],basis=['e'],dtype='scalar')
  ```
  Note that chosing the scalar type the defining a multivector basis other than the scalar basis, returns an error. For example the following will return an error
  ```python
  x = ga.multivector([[1.234],[2.345]],basis=['e1'],dtype='scalar')
  ```

### Using numpy to create multivectors and multivector arrays

```python
import gasparse
import numpy as np
n = 10
ga = gasparse.GA(3)
arr = np.random.rand(n,ga.size(1,2)) # innermost dimension must be the the size of grades 1 and 2
x = ga.multivector(arr.tolist(),grades=[1,2]) # only accepts lists as input
```

### Coding Strategy
1. Set errror strings only in the outermost call (the first function that is called from python)
Still trying to averiguate the best strategy for error setting...

### Usefull Commands
- Leak check:
  `valgrind --leak-check=full --tool=memcheck --suppressions=valgrind-python.supp python3 -E snippets/leak_check_script.py`
- Debug python code with gdb
  `gdb -ex r --args python3 snippets/test_template_gen.py`

### Generating Code and Compiling
Before compiling the C code ensure that some packages are installed

`sudo apt-get install python3-dev`
```shell
cd sparse-multivectors
python3 code_gen_large.py # Run if a change in multivector_large.c.src is made
python3 genalgebra.py # Run if added a new algebra or changed multivector_gen.c.src
python3 setup.py build # builds gasparse
python3 setup.py build --genalgebras # generates the geometric algebras operations
```

Compiling and building the package to install with pip
```shell
cd sparse-multivectors
python3 setup.py build
python3 setup.py sdist bdist_wheel
pip install . # Installs the contents of the current directory
```


### Code Structure

`gasparse.c` -> algebra declarations and initializations

`multivector_object.c` -> multivector array object calls

`multivector_types.c, multivector_large.c, multivector_gen.c` -> Different implementations of the different computation modes. 


### Note
When creating specific function types if we want that the compiler warns incompantible function pointers do not cast a function pointer otherwise if the function is not as specified in the typedef we will probably have segmentation faults.

## BIG TODO

- `DenseMultivector` operations using position and sign tables.
- `GradedDenseMultivector` operations using position and sign tables, for each grade of the multivector 
  + binary `table[grade_left][grade_right][grade_out]`
  + ternary ``table[grade_left][grade_center][grade_right][grade_out]``
  + Only create tables when needed
  + Use a flyweight pattern to decrease memory consumption
- Faster Scalar products `(a*b)(0)` -> `mv.sp(a,b)` `mv.sp` dispatches the scalar product
- Change the name of `BladesMultivector` to `GradedSparseMultivector`
- Implement the graph/linked list strategy for the `BladesMultivector` in the `"large"` computation mode

## SMALL TODO
1. `ga.multivector`: 
	- $\checkmark$ Check the size of the innermost dimension. It must have the same size as the specified basis or grades.
	- $\checkmark$ Check if the specified grades are valid.
	- scalar types do not need to specify grades or basis.
1. Grade projection also as `a(1,2,3)` and `a([1,2,3])` (read the tuple directly) $\checkmark$
1. Get sizes with `ga.size(1,2)` and `ga.size([1,2])` $\checkmark$
1. Change code for `x.list()` and `ga.size()` making it accept variable number of arguments like `x.list(1,2)` instead of a list like `x.list([1,2])` $\checkmark$
1. Indexing multivector arrays.
    1. $\dagger$ single index `x[8][9]`  $\checkmark$
    1. Multiple indices `x[8,9,10]` $\checkmark$
    1. A simple slice `x[:,1]` $\checkmark$
    1. Range of indices `x[10:15]`
    1. $\dagger$ Check if the indices are in range of the array
1. Masking multivector arrays. I want to be able to eliminate some multivectors that do not satisfy some requirements.
1. Grade Einsum
1. The size of dense type arrays do not need to be inside of the struct object (Extract it from the algebra).
1. Grade projection should accept the grade bool array instead of the grades 
1. $\dagger$ $\checkmark$ Add option to print multivectors with `*` instead of `^` to be able to copy it's value and declare it as variable in python, otherwise because of the operator precedence the output would not result in what is expected, that is these are not equivalent
    ```python
    x = 1.234*e1 + 5.342*e4 + 5.678*e12
    x = 1.234^e1 + 5.342^e4 + 5.678^e12
    ```
    to be equivalent we would need to put parenthesis around each operation as (or just add parenthesis)
    ```python
    x = (1.234^e1) + (5.342^e4) + (5.678^e12)
    ```
1. $\dagger$ When printing big multivector arrays make sure that the characters do not surpass the maximum size. (Use the approach used by numpy (the dots))
1. Write a script to test all of the functionalities of the package (Also include valgrind to check for memory leaks) The common cases where we usually have leaks is when there is some error, check for leaks in such cases. Use
  ```python
  try:
    # lines of code
  except Exception:
    pass
  ``` 
  to ignore exception errors.
1. $\dagger$ Ignore multivector elements via relative difference, instead of just small absolute value (make that as an option for the user, consider elements of the same grade) $\checkmark$
4. $\dagger$ Change the value for which multivectors get converted to zero. Should be an option when generating the algebra (epsilon) $\checkmark$
8. $\dagger$ Give an option to not remove relatively small multivector elements. Set the precision to zero. $\checkmark$
1. Make `x.grade()` retrieve the grades based on the relative difference between the basis elements of the different grades 
1. Implement cast to some algebra (Add the option to project to that algebra)
10. Implement grade involution (Important for reflection)
11. Implement integer powers of multivectors
12. Print multivectors by grade
13. Generate random multivector arrays
20. Implement einsum where the product is some specified product (see [general product]("general product"))
22. Implement the above in **CUDA**
23. Create a python extension from the **C** and **CUDA** code
25. Compile for diferent versions of python
1. $\dagger$ Documentation

- Pretty printing multivector arrays:
  + $\checkmark$ Add space before square brackets (align brackets that are of the same dimension) 
    ```python
    [[1e-32*e1 + 0*e2 + 0*e3],
     [0*e1 + 2*e2 + 0*e3],
     [0*e1 + 0*e2 + 3*e3]]
    ```
  + Better printing for multiple dimensions $\checkmark$
  + Limit the amount of multivectors printed (use `...` to print big arrays)
  + Printing multivector arrays using realloc. $\checkmark$
  + Option to print with metric. $\checkmark$
  + Single multivectors not in quare brackets `[[1*e1],[2*e1]]` $\rightarrow$ `[1*e1,2*e1]`

1. $\checkmark$ Check indices when iterating mvarrays `x[k]`
1. $\checkmark$ The following does not output correctly
```python
arr = np.random.rand(3,3,ga.size(1,2)) 
x = ga.mvarray(arr.tolist(),grades=[1,2])
```
but printing `x[0]`, `x[1]` and `x[2]` yields correct results.


1. Mapping protocol
	- $\checkmark$ Indexing multivector arrays `x[i]`
	- $\checkmark$ Getting length and shape from multivector arrays `len(x)` and `x.shape()`
	- setting values to items at some index `x[5] = y`
  - indexing with slices `x[2:]`

1. $\checkmark$ Solve segmentation fault when asking for size `ga.size(1)`

1. the `unary_sparse_scalaradd` involves looking at the precision, need to review this function.


<!-- 1. Scalar type multivectors are agnostic to the grade of the multivectors as such returning  -->

**$\dagger$ After these are done make the repo public!!**

### Other operations:

Except logarithm and exponential these operation can be computed using the general product

+ **Grades Projection** - $\langle y\rangle_{j_1,j_2,\dots}$
+ **Inner Product**  - $x \cdot y$
+ **Wedge Product** - $x \wedge y$
+ **Reverse** - $x^\dagger$
+ **Sandwich Product** - $xyx^{-1}$
+ **Projection** - $x\cdot y y^{-1}$
+ **Rejection** - $x\wedge y y^{-1}$
+ **Exponential** - $\exp(x)$
+ **Logarithm** - $\text{Log}(x)$
+ **Inverse** - $x^{-1}$
+ **Norm** - $\langle xx^\dagger\rangle$
+ **Dual** - $x^* = xI^{-1}$
+ **Meet** - $m = x\cap y\equiv (x\cdot m^{-1})\wedge y$
+ **Join** - $j = x\cup y\equiv (x\cdot j^{-1})\cdot y$



## 1. General Product

grade selection $\rightarrow$ geometric product $\rightarrow$ grade selection.

$\langle\langle x\rangle_{i_1,i_2,\dots}\langle y\rangle_{j_1,j_2,\dots}\rangle_{k_1,k_2,\dots}$

where $\langle \cdot\rangle_{i_1,i_2,\dots}$ selects the grades $i_1,i_2,\dots$

## 2. Data-types

### 2.1 Sparse

- Array of values with associated bitmap.
- Each bitmap represents each basis blade.

```c
typedef struct sparse{
    int *bitmap;
    float *value;
    unsigned int size;
}sparse;
```

### 2.2 Grade-Sparse

- Sparse multivectors organized by grade.
- Each grade element has a sparse representation.

```c
typedef struct blades{
    sparse *data;
    unsigned int *grade;
    size_t size;
}blades;
```

### 2.3 Dense

- Has the full representation of the multivector the index is its own bitmap.
- The size is the size of the algebra $2^{p+q+r}$.

```c
typedef struct dense{
    float *value;
    unsigned int size;
}dense;
```




### 2.4 Pseudo Versors

**Vector Types**
```c
typedef struct dense_vector{
    float *value;
}dense_vector;

typedef struct sparse_vectors{
    float *value;
    unsigned int size;
}dense_vectors;
```

Data structures for multivectors that are expressed as the product of vectors:
```c
typedef struct product_multivector{
  vector *vectors;
  unsigned int size; // How many vectors
  ProductType ptype; // The type of product (usually geometric product and outer product)
}
```

Blades can be allways written as $A = a_1\wedge a_2\wedge \cdots \wedge a_r$ and if the $a_i$ are orthogonal then we can drop the wedge and simply write $A = a_1 a_2 \cdots  a_r$. Operations under this representation are harder to employ but some simplifications might appear and we might also have some numerical benefits. Blades are guaranteed to be blades, the space that it takes to store can be smaller then if using the geometric algebra basis. Versors $A$ can allways be written as the geometric product of not nececerily orthogonal vectors $a_i$ as $A = a_1 a_2 \cdots  a_r$.

Using this representation we for each sandwich with a vector $axa^{-1} = 2x\cdot aa^{-1} - x$ yields $4n$ operations $n$ for the inner product $n$ more for the product of the scalar $x\cdot a^{-1}$ with the vector and $2n$ more for the subtraction.  

```c
typedef struct versor{
    vector *vectors; // Array of vectors (sparse or dense vectors)
    void *multivector; // A general multivector
    unsigned int size;
    VersorSubType Type; // SubType of versor: blade, ortho blade, prod gen vectors, multivector
    VersorDataRepr Repr; // Type of representation: sparse, blades (only for multivector), dense
}versor;
```

This type of representation is used to express linear transformations $f:\mathcal{A}_{p,q}\mapsto \mathcal{A}_{p,q}$
Multivectors for which $AA^\dagger$ is a scalar and that $AxA^\dagger$ is grade preserving.

Subtypes:
- Outer product of vectors $A = a_1\wedge a_2\wedge\cdots\wedge a_r$
- Product of orthogonal vectors $a_ia_j = a_i\wedge a_j$, for $j\neq i$. $A = a_1\wedge a_2\wedge\cdots\wedge a_r$
- Product of general vectors $A = a_1\wedge a_2\wedge\cdots\wedge a_r$. Where $a_ia_j\neq 0$ for $i\neq j$.
- Multivector $A: AA^\dagger = \langle AA^\dagger\rangle$ 

In all of the above subtypes we can have that $AA^\dagger=0$. That is, $A$ can be a null vector.

## Representation of Linear Transformations
Let $f:\mathcal{A}_{p,q}\mapsto \mathcal{A}_{p,q}$ and let $x\in\mathcal{A}_{p,q}$

- Matrix representation $f(x) = \sum_{ij}f_{ij} x\cdot e^ie_j $
- Multivector representation: $f(x) = \langle F(x)\rangle_1$
  1. Geometric Product: $F(x) = \sum_{k} A_kxB_k$
  1. Inner Product: $F(x) = \sum_{k} x\cdot A_kB_k$
  1. Outer Product: $F(x) = \sum_{k} x\wedge A_kB_k$

When $B_k = \alpha_k A_k^\dagger$ simplifications arise: Specifically When the multivectors $A_k$ are certain types of multivectors then we can drop the grade projection, particularly if $A_k$ is a blade in 1.,  2. or 3. or when the $A_k$ are Pseudo Versors in 1. we have $F(x) =  \langle F(x)\rangle_1$. Note that representations 2. and 3. can be expressed with 1. by relating inner and outer products with geometric products. Note that for 2. when $A_k$ and $B_k$ are vectors this is the usual vector representation of a transformation $f(x) = \sum_{k} x\cdot a_kb_k$.

### Computing Transposes of Linear and Multilinear Transformations

$\bar{F}(X) = \partial_{Y}^\dagger X^\dagger*F(Y)$ implies that $F(X)*Y^\dagger = X^\dagger * \bar{F}(Y)$

### Computing inverses of Linear Functions
$f^{-1}(A') = \bar{f}(A'I')[\bar{f}(I')]^{-1}$

### Computing Outermorphisms of linear Functions
$\underline{f}(X_r) = X_r*\partial_{(r)}f_{(r)}$

$\bar{f}(X_r) = \partial_{(r)}X_r*f_{(r)}$



**Implementing the Linear and Multilinear Operators/Functions**

Operations on vectors can be simplified greatly. Specifically the product $\langle AxB\rangle_1$ can be computed as a ternary operation.

## Linear and Multilinear Functions/Operators
Defining functions(MultiLinear Operators):
```python
F = gasparse.multilinear.operator(A,B,grades=3) # Create a grade preserving multilinear operator. F(X) = sum_i (A[i]*X(3)*B[i])(3)
F = gasparse.multilinear.operator(F_array,grades=(1,2)) # Create a multilinear operator for grades 1 and 2 from an array F(X) = sum_ij F_ij (X|fi)*fj
F = gasparse.multilinear.operator(A,B) # Create a general multilinear operator F(X) = sum_i A[i]*X*B[i]
F = gasparse.multilinear.operator(A,B,linear=True) # An operator from vectors to vectors: F(X) = sum_i (A[i]*X(1)*B[i])(1)
F = gasparse.multilinear.operator(A,right=True) # An operator that applies a binary operation
F = gasparse.multilinear.operator(A,right=True,grades=1) # F(X) = (X(1)*A)(1) usefull to define skew transformations
```


## Grade Preserving Operations and outermorphisms of scaled orthogonal Functions

When $f(x) \equiv \mu UxU^\dagger$ with $UU^\dagger=1$ then the outermorphism of $f$ takes the form: $\underline{f}(X_r) = \mu^r UX_rU^\dagger$

$\boxed{F(X) = \langle A\langle X\rangle_r B\rangle_r}$

Calls to specific operations: (Generate the ternary operation tables)

- $aX_ra$
- $axa$
- $AX_rA$
- $a_1a_2\cdots a_rXa_ra_{r-1}\cdots a_1$
- $AXA$ where $A$ is a versor


The most generic ternary operation on $X$ is $\langle A\langle X\rangle_{r_1,r_2,...,r_\ell} B\rangle_{s_1,s_2,...,s_t}$

**TODO (MultiLinear Operators)**
All the operations done inside of the operator must be hidden from the end user.

- I want to be able to apply an operator to a multivector array.
- I want to be able to define multiple operators simultaneously from a multidimensional array of multivectors.
- The multilinear operator can have a decomposition `gasparse.multilinear.eig(F)` (the eigendecomposition of the operator F, returns eigenvalues and eigenblades) 
- Decompositions for normal and orthogonal transformations
- Linear Operators as a binary operation `F = gasparse.multilinear.operator(A,right=True)` $F(X) = XA$. This is usefull to define skew transformations $F(x) = x\cdot B$ where $B$ is a bivector.
1. Linear Operators -- grade one preserving Operators: Generate the ternary operator table for $axa$, $AxA^\dagger$ and $\langle AxA^\dagger \rangle_1$
1. Overload calls to the operator Object `F(X)` applies the linear operator `F` to the multivector `X`.



## Questions

**Hashmap/Hashfunction** Using an hashmap/hashfunction to convert bitmaps to grade, position.
What is an efficient way to compute hashmaps/hashfunctions?

What is the outermorphism of $f(x) = UxU^\dagger$ when $UU^\dagger=0$. Is it zero? $f(x)\wedge f(y) = (UxU^\dagger)\wedge(UyU^\dagger) = \tfrac{1}{2} \left( UxU^\dagger UyU^\dagger -  UyU^\dagger UxU^\dagger \right) = 0 $

**Cleaning Sparse Multivectors**
For the sparse representation of multivectors should I have code that the user (python code) executes to clean sparse multivectors??
How should I 'clean' the sparse multivectors. Should I compare the greatest value with the all the other values. distance between values `abs(a[i]-a_max) < eps` or `abs(a[i]/a_max) < eps` or compare with some percentage of max `abs(a[i]) < abs(a_max)*eps`. `eps` is a value smaller then one.

What is the computation complexity difference between computing $UXU^\dagger$ and $u_1u_2\cdots u_m X u_m u_{m-1}\cdots u_1$?? 

### Ternary Unique Grade Operations
 - Create a new field for the structs holding the operations for the different types
 - Store each multivector with its own grade.
 - Specialized products: 
    - (Unique,Unique,Unique) -> Unique
    - (General,Unique,General) -> Unique
    - (General,General,General) -> General
 - (There is a cost associated with this type of operation) I need to convert bitmaps to position and grade
 - In the future create a new type of multivector which is indexed by grades and positions, find a efficient method to compute positions when multiplying mvs in that data format.
 - For generated types create an array of functions for each grade: `func[grade0][grade1][grade2][grade_out]` then iterate over all the grades of the multivector (generated types are already indexed by grade and position) Note that for 3D CGA this amounts to $6^4 = 1296$ functions for the ternary product. Also a lot of this functions are going to be empty. Need a way to indicate that a grade is to be ignored...

***Question: How can optimize the sandwich product with versors?*** 

Versors are either even or odd multivectors. They do not have to necessarily span all even or all odd grades. The optimization can be done for the multivector by which the multilinear operation $F(X) = VXV^\dagger$ is applied to. If $X$ is of unique grade then we can surelly define a table for all even grade and for the grade of $X$. Table for grade preserving multilinear transformation: Define a grade wise ternary map `map[grade][even/odd]` where `grade` is the grade of the output and of the input multivector. The table can then be indexed by grade/position or by bitmaps.


A specialized function for the sandwich products $\alpha AXA^\dagger$, $\alpha X\wedge AA^\dagger$, $\alpha X\cdot AA^\dagger$.

The function `specialized_ternary_blades_product0` is a specialized ternary function that computes the products $\alpha A\odot B\odot C$ or $\alpha A\odot B\odot C^\dagger$ where $\odot$ means any of the different products geometric, inner, outer, regressive, etc..

Creating a subpackage for multilinear operations:
 - Ternary operations.
 - Creating functions from arrays of multivectors.
 - Functions applied to multivectors execute the ternary operations.
 - Need to have:
    1. Ternary maps for the different products.


How do I add it as an optional subpackage?? 
Since ternary maps occupy a lot of space can I make them be only computed when necessary?
Creating the first in code multilinear function should automatically compute the ternary table. 
Or after executing the ternary product checks if the table is already built.

**Only generate tables that are needed for that operation!!**

How to generate tables at execution time for a specific product?
How expensive is that operation?

Note that memory for the ternary product is quite expensive since `map[ptype0][ptype1][grade0][grade1][grade2][gradeout]` in 3D CGA has $4^2*6^4=20736$ entries, which by itself is quite expensive.

### Code generating ternary products
The type of multivector that will benifit the ternary products are the code generated types. While generating a lot of ternary functions might create a very long `c` file the execution time might drop dramatically.

### A less memory expensive way to deal with ternary products
- Use the binary sign maps.
- Use the bitmaps to compute the resulting basis blades.


## Multilinear Functions
### TODO
  1. Implement a fast matrix computation algorithm using binary products.
      1. Generated types (indexed by grade,position)
      1. Dynamically allocated types (indexed by bitmap)
  1. Implement new type indexed by grade and position 
      1. Implement the graded products (binary and ternary)
      1. Generate code for ternary and binary graded products   
  1. Implement the eigendecomposition for multiplicity greater then one.


### Future Optimizations

- The scalar product (Run through the diagonal elements of the cayley table) efficient for dense and generated multivectors
- `C` compiler arguments.
- copying multivector arrays using memcopy
- `'large'` multivector array type:
  + signs in a trenary representation, each triac of a `long` or a `float` represents the sign
  + Compute product while maintaning a sparse output. That is do not alloc memory for the entire multivector output. 
  Use a graph or something like that and do a search. 
- Using `memset` and `calloc` to atribute initial values to the values and bitmaps 
  
  ```c
    dense->bitmap = calloc(size,sizeof(Py_ssize_t));
    dense->value = calloc(size,sizeof(ga_float));
    memset(dense->bitmap, -1, size * sizeof(Py_ssize_t));
  ```



To avoid iterating over the entire array of multivectors:
  - Store the indices in the smallest data structure possible 
  - Remove repeated indices
  - Atribute the value and bitmap for each index