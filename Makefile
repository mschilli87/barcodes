#######################
# general information #
#######################

# file:       Makefile
# created:    2016-01-28
# author(s):  Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose:    automate building of C program for barcode collapsing for drop-seq data


######################################
# change log (reverse chronological) #
######################################

# 2016-01-28: initial version (thinned out version of corresponding project Makefile)


#########################
# parameter definitions #
#########################

# command line flags to use for gcc
GCC_FLAGS:=-Wall -Ofast -Wpedantic -Werror -march=native -mtune=native


####################
# path definitions #
####################

# (absolute) path of this Makefile
MAKEFILE:=$(realpath $(lastword $(MAKEFILE_LIST)))

# (absolute) path of the directory containing this Makefile
MAKEFILE_DIRECTORY:=$(dir $(MAKEFILE))

# (absolute) path of C source file for program to re-map cell barcodes
REMAP_BARCODES_C:=$(MAKEFILE_DIRECTORY)remap_barcodes.c

# (absolute) path of binary file of program to re-map cell barcodes
REMAP_BARCODES_BINARY:=$(REMAP_BARCODES_C:.c=)


################
# make options #
################

.DELETE_ON_ERROR :
.SUFFIXES:
.SECONDARY:


#######################
# command definitions #
#######################

# command used to run gcc
GCC:=gcc

# command used to compile C source files to target executable file
COMPILE_FROM_C_SOURCE:=$(GCC) $(GCC_FLAGS) -o


##################
# common targets #
##################

# if no target was specified, build program to remap barcodes
.PHONY : all
all : $(REMAP_BARCODES_BINARY)


#################
# build program #
#################

# compile program to remap barcodes from C source file
$(REMAP_BARCODES_BINARY) : $(REMAP_BARCODES_C) | $(dir $(REMAP_BARCODES_BINARY))
	$(COMPILE_FROM_C_SOURCE) '$@' '$<'
