# Introduction

KD-Tree is a standard data structure for indexing data, especially in 3D space. It is an extension of Binary-Space Partition (BSP) to more than one dimension. For more information on KD-Tree, please refer Wiki.

This repository is a Fortran implementation of KD-Tree. We need KD-Tree in scientific HPC scenarios, so Fortran is still the right language, and also we need more modern library interfaces.

This repository is a fork from [Li Dong's one](https://github.com/dongli/fortran-kdtree). Trees are built with real(4) data, as opposed to the initial repository.

# Compilation

To get and compile the library,

1. Clone the repository:

```bash 
git clone https://github.com/dvallesp/fortran-kdtree
```

2. Move to the repository folder and compile the library:

```bash
cd fortran-kdtree
cmake CMakeLists.txt
make
```

3. Move the `libfortran_kdtree.a` and `kdree.mod` file to the folder where you want to use the library.

4. Compile your code `yourcode.f` including the library:

```bash 
    gfortran -O3 -fomit-frame-pointer (-whatever-option-you-might-add) yourcode.f libfortran_kdtree.a -o your_program.x
```

# Usage

```fortran
  use kdtree

  real(4), allocatable :: x(:,:)       ! num_dim, num_point
  integer, allocatable :: ngb_idx(:) ! num_ngb
  type(kdtree_type) tree
  integer i,j
  real(4) a(num_dim)

c   BUILDING A TREE
  ! Some random data
  do i=1,num_point
  do j=1,num_dim
    call random_number(a)
    arr(j,i)=a
  end do
  end do

  ! Tree built
  call tree%build(x)

c   SEARCH IN THE TREE
  ! Some random point
  call random_number(a)
  num_neigh=32
  allocate(ngb_idx(num_neigh))
  call tree%search(a, ngb_idx)
  ! Now ngb_idx contains the indices of the num_neigh first neighbours.
```

Note that the calls to `tree%search()` can be in principle OMP parallelised. In my tests, there is an almost-linear performance improvement from 1 to 4 threads, but not much beyond it.

# Contributors

- Li Dong ([original repository](https://github.com/dongli/fortran-kdtree))
- David Vallés-Pérez (this fork, which just changes building data to `real(4)` and adds a bit of information).
