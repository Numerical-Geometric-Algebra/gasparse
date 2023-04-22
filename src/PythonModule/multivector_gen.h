#ifndef MULTIVECTOR_GEN_H_
#define MULTIVECTOR_GEN_H_

#define N_GEN_SUBTYPES 2

// enum with all the different types
typedef enum {
  gen_MultivectorTypeMIN = 2,
  gen_MultivectorType_dense0 = 3,
  gen_MultivectorType_blades0 = 4,
  gen_MultivectorTypeMAX} gen_MultivectorType;

// code generate the size of this array
extern PyMultivectorSubType gen_subtypes_array[N_GEN_SUBTYPES];

#endif // MULTIVECTOR_GEN_H_