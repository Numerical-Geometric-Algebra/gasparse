# sparse-multivectors

Creating a python library to do computations using sparse representation of multivectors.

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

Note that the last two lines generate the same algebra yet the symbols represent different basis, as such they are not compatible. Furthermore the last line is not compatible with the vga algebra, since the first basis element is negative. Note that the following will result in an **error**:
```python
vga = gasparse.GA(3)
cga = gasparse.GA(metric=[-1,1,1,1,1])
t = cga.multivector([1,2,3,4,5],grades=1)
y = vga.multivector([1,2,3],grades=1)
z = t + y
```

- Initialize a basis
```python
basis = cga.basis()
locals().update(basis) # Update all of the basis blades variables into the environment
```
The basis will be available as `e1`, `e2`,..., `e12`,...,`e12345`. 

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
  - For multiple grades `ga.size([1,2])`
  - For the whole algebra (empty argument) `ga.size()`

- Zero grade projection now outputs objects of type float instead of `gasparse.multivector`
- The user can now use `x.grade()` to get the grade of a multivector if it has multiple grades then it returns `None`

- The function list accepts grades as input, can be a single integer or a list 
```python
values,basis = x.list([0,1]) # returns grades zero and one
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
```
Note that the division `/` operator is not defined since for general multivectors is not simple to compute. Instead for homogeneous multivectors we consider the operation
```python
(1/((x*~x)(0)))*x
```
(Maybe add the option to be able to divide by scalars...)

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
 x = ga.multivector(list)
 print(type(x(0)))
 ```
should be `float`.

### Coding Strategy
1. Set errror strings only in the outermost call (the first function that is called from python)
Still trying to averiguate the best strategy for error setting...

### Usefull Commands
- Leak check:
  `valgrind --leak-check=full --tool=memcheck --suppressions=valgrind-python.supp python3 -E snippets/leak_check_script.py`
- Debug python code with gdb
  `gdb -ex r --args python3 snippets/test_template_gen.py`

- Generate code and compile
```shell
cd sparse-multivectors
python3 code_gen.py # Run if a change in multivector.c.src is made
python3 code_gen_large.py # Run if a change in multivector_large.c.src is made
python3 genalgebra.py # Run if added a new algebra or changed multivector_gen.c.src
python3 setup.py build # builds gasparse
python3 setup.py build --genalgebras # generates the algebra
```

## TODO
1. For operations involving `PyMultivectorObject` that are transparent to the subtypes, it might make sense to create a new empty multivector and then pass that as reference to function that operates directly on the type.
1. Test mixed type operations involving the generated code
1. Test the cast function
1. $\checkmark$ Scalars should output as floats and not as gasparse.multivector objects
1. $\checkmark$ Function to check grade (return `-1` or `None` if it isnt of unique grade)
1. Write a script to test all of the functionalities of the package (Also include valgrind to check for memory leaks) The common cases where we usually have leaks is when there is some error, check for leaks in such cases. Use
```python
try:
   # lines of code
except Exception:
   pass
```
to ignore exception errors.
1. Ignore multivector elements via relative difference, instead of just small absolute value (make that an option for the user, consider elements of the same grade)
1. $\checkmark$ Read and Write Multivectors by grade
   
   - [ ] Read/write without defining basis (just defining grade)
   - [ ] Read/write in some basis (can be given as string or as an arbitrary basis in some algebra)
     - [ ] Might need to compute the reciprocal basis (or give that basis as input)
   - [ ] Option to also output the basis as bitmap or as the basis of those grades of that algebra
     - [ ] The basis can be multivectors, or positive integers (bitmaps)
3. Return multivectors in a list by a specified grade
4. Change the value for which multivectors get converted to zero. Should be an option when generating the algebra (epsilon)
1. Change code for `x.list()` and `ga.size()` making it accept variable number of arguments like `x.list(1,2)` instead of a list like `x.list([1,2])`
1. Implement cast to some algebra (Add the option to project to that algebra)
6. Numpy integration
   
   - [ ] Read numpy arrays
   - [ ] Generate random multivectors by grade or by given basis
   - [ ] Numpy Array to array of multivectors (by a given basis, by grade or by all elements of the algebra)
   - [ ] Output multivectors as a numpy array by grade (also output basis, consider a basis as input)
7. Do the same but also with lists of lists
8. Give the option to compile the code without the zero value checking (checking if it is close to zero)
9. Write code to retrieve multivectors as lists
   
   - [ ] Multivectors by grade
   - [ ] Also output bitmap
   - [ ] Convert data to and from numpy
   - [ ] Convert multilinear and linear transformations into matrices
10. Implement grade involution (Important for reflection)
11. Implement integer powers of multivectors
12. Print multivectors by grade
13. Generate multivector
    
    - [ ] Random by grade
    - [ ] From any list
    - [ ] Multiple grade selection from multivectors
14. Multivector arrays!!!!!!
15. Implement the geometric product in all data representations in **C**
    
    - [ ] Data type conversions
      + [x] sparse to grade sparse
      + [x] sparse to dense
    - [x] geometric product sparse
    - [x] geometric product grade-sparse
    - [x] geometric product dense
      + [x] Direct mapping
      + [x] Inverse mapping
16. Implement grades selection for each data type $\langle y\rangle_{j_1,j_2,\dots}$
    
    + [x] dense
    + [x] sparse
    + [x] grade-sparse
17. $\checkmark$ Implement the general product
18. $\checkmark$ Implement sum for each type (still need to test):
    
    + [x] append
    + [x] add
19. $\checkmark$ Generate tables for the outer and inner products
20. Implement einsum where the product is some specified product (see [general product]("general product"))
21. Use jinja to generate code for different types of values
22. Implement the above in **CUDA**
23. Create a python extension from the **C** and **CUDA** code
24. Integrate this library into pytorch
25. Compilar para diferentes versoes de python

### Other operations:

Except logarithm and exponential these operation can be computed using the general product

+ **Grades Projection** - $\langle y\rangle_{j_1,j_2,\dots}$
+ **Inner Product**  - $x \cdot y$
+ **Wedge Product** - $x \wedge y$
+ **Sandwich Product** - $xyx^{-1}$
+ **Projection** - $x\cdot y y^{-1}$
+ **Rejection** - $x\wedge y y^{-1}$
+ **Exponential** - $\exp(x)$
+ **Logarithm** - $\text{Log}(x)$
+ **Inverse** - $x^{-1}$
+ **Norm** - $\langle xx^\dagger\rangle$
+ **Reverse** - $x^\dagger$
+ **Dual** - $x^* = xI^{-1}$
+ **Meet** - $m = x\cap y\equiv (x\cdot m^{-1})\wedge y$
+ **Join** - $j = x\cup y\equiv (x\cdot j^{-1})\cdot y$

Before compiling the C code ensure that some packages are installed

`sudo apt-get install python3-dev`

### Code

- [ ] change all `unsigned int size` to `size_t size`

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
}sparse;
```

