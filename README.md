# sparse-multivectors

Creating a python library to do computations using sparse representation of multivectors.

## TODO

1. Read and Write Multivectors by grade
   
   - [ ] Read/write without defining basis (just defining grade)
   - [ ] Read/write in some basis (can be given as string or as an arbitrary basis in some algebra)
     - [ ] Might need to compute the reciprocal basis (or give that basis as input)
   - [ ] Option to also output the basis as bitmap or as the basis of those grades of that algebra
     - [ ] The basis can be multivectors, or positive integers (bitmaps)
2. Scalars should output as floats and not as gasparse.multivector objects
3. Change the value for which multivectors get converted to zero. Should be an option when generating the algebra (epsilon)
4. Function to check grade (return `-1` or `None` if it isnt of unique grade)
5. Numpy integration
   
   - [ ] Read numpy arrays
   - [ ] Generate random multivectors by grade or by given basis
   - [ ] Numpy Array to array of multivectors (by a given basis, by grade or by all elements of the algebra)
   - [ ] Output multivectors as a numpy array by grade (also output basis, consider a basis as input)
6. Do the same but also with lists of lists
7. Give the option to compile the code without the zero value checking (checking if it is close to zero)
8. Write code to retrieve multivectors as lists
   
   - [ ] Multivectors by grade
   - [ ] Also output bitmap
   - [ ] Convert data to and from numpy
   - [ ] Convert multilinear and linear transformations into matrices
9. Implement grade involution (Important for reflection)
10. Implement integer powers of multivectors
11. Print multivectors by grade
12. Generate multivector
    
    - [ ] Random by grade
    - [ ] From any list
    - [ ] Multiple grade selection from multivectors
13. Multivector arrays!!!!!!
14. Implement the geometric product in all data representations in **C**
    
    - [ ] Data type conversions
      + [x] sparse to grade sparse
      + [x] sparse to dense
    - [x] geometric product sparse
    - [x] geometric product grade-sparse
    - [x] geometric product dense
      + [x] Direct mapping
      + [x] Inverse mapping
15. Implement grades selection for each data type $\langle y\rangle_{j_1,j_2,\dots}$
    
    + [x] dense
    + [x] sparse
    + [x] grade-sparse
16. $\checkmark$ Implement the general product
17. $\checkmark$ Implement sum for each type (still need to test):
    
    + [x] append
    + [x] add
18. $\checkmark$ Generate tables for the outer and inner products
19. Implement einsum where the product is some specified product (see [general product]("general product"))
20. Use jinja to generate code for different types of values
21. Implement the above in **CUDA**
22. Create a python extension from the **C** and **CUDA** code
23. Integrate this library into pytorch

1. Compilar para diferentes versoes de python

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

