#ifort bemrtm.f  bipts.f  discrete.f  parameter.f  -o bemrtm
FC=gfortran
FLAGS=

all: bemlayer 

bemlayer: bemlayer.f  bipts.f  discrete.f  fmat1.f  fmat2.f  fmats.f  FRHM.f  gasave.f  ludcmp.f  parameter.f
	$(FC) $(FLAGS) bemlayer.f  bipts.f  discrete.f  fmat1.f  fmat2.f  fmats.f  FRHM.f  gasave.f  ludcmp.f  parameter.f -o bemlayer


