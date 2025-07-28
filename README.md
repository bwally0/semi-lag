# Semi-Lagrangian 3rd Order Model

- Mohamed Moustaoui, PhD
- Matei Georgescu, PhD

## intel compiler

To compile the code using the Sol intel compiler
you need to load the intel compiler:

``` bash
module load intel/intel-oneapi-2022.1.0
```

To compile use the command:
``` bash
ifort -FR -o semi_lag semi_lag_3th_order_free.f
```

To execute use the command:
``` bash
./semi_lag
```

The code run for 3 days
The output for day 3 is written to a file with a frequency of 1 hour
The name of the file generated is: `temp_fdx.dat`
It contains the 2D fields of temperature for each hour.

### post processing

For postprocessing compile and execute the codse `tmin_tmax.f`:
``` bash
ifort -FR -o tmin_tmax tmin_tmax.f
./tmin_tmax
```

The files `tmin_max_dt_tav_fdx.dat` and `phx_temp_fdx.dat` will be generated.

1. `tmin_max_dt_tav_fdx.dat` is a binary file that contains the 2D fields
of tmin, tmax DT (diurnal range) and tav (average temperature)

2. `phx_temp_fdx.dat` is an ascii file that you can open and read.
It containes 2 colomns: hour of the day  and temperature at PHX airport.

The expected outcome should be the same as the file `phx_temp_fdx.dat_save`
which is located in the same directory

## gnu compiler

If you want to use the gnu compiler edit the files:
`semi_lag_3th_order_free.f` and `tmin_tmax.f` Dand do the following changes:

look for the lines:

``` Fortran
   integer, parameter :: recl_f=1*nx*ny  ! for intel ifort
!   integer, parameter :: recl_f=4*nx*ny  ! for gnu gfortran
```

and replace the with:

``` Fortran
!   integer, parameter :: recl_f=1*nx*ny  ! for intel ifort
   integer, parameter :: recl_f=4*nx*ny  ! for gnu gfortran
```

To compile and execute `semi_lag_3th_order_free.f` use the commands:
``` bash
gfortran -ffree-form -o semi_lag semi_lag_3th_order_free.f
./semi_lag
```

To compile and execute `tmin_tmax.f` use the commands:
``` bash
gfortran -ffree-form -o tmin_tmax tmin_tmax.f
./tmin_tmax
```
