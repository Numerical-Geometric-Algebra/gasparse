# sparse-multivectors

Creating a python library to do computations using sparse representation of multivectors.

## DONE

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
  Note that the values can be either integer or float. (Everything gets cast to `ga_float`.)
- Function that returns sizes:
  - For each grade `ga.size(2)`
  - For multiple grades `ga.size([1,2])`
  - For the whole algebra (empty argument) `ga.size()`

- Zero grade projection now outputs objects of type float instead of `gasparse.multivector`
- The user can now use `x.grade()` to get the grade of a multivector if it has multiple grades then it returns `None`

- The function list accepts grades as input, can be a single integer or a list 
```python
values,bitmaps = x.list([0,1]) # returns grades zero and one
values,bitmaps = x.list() # returns all elements
```
The first variable of the output are the values atributed to each basis, the second are the bitmaps that indicates the basis (Maybe it would be better if it where a basis in a multivector type), 

### Coding Strategy
1. Set errror strings only in the outermost call (the first function that is called from python)

## TODO

1. $\checkmark$ Scalars should output as floats and not as gasparse.multivector objects
1. $\checkmark$ Function to check grade (return `-1` or `None` if it isnt of unique grade)

1. $\checkmark$ Read and Write Multivectors by grade
   
   - [ ] Read/write without defining basis (just defining grade)
   - [ ] Read/write in some basis (can be given as string or as an arbitrary basis in some algebra)
     - [ ] Might need to compute the reciprocal basis (or give that basis as input)
   - [ ] Option to also output the basis as bitmap or as the basis of those grades of that algebra
     - [ ] The basis can be multivectors, or positive integers (bitmaps)
3. Return multivectors in a list by a specified grade
4. Change the value for which multivectors get converted to zero. Should be an option when generating the algebra (epsilon)
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

