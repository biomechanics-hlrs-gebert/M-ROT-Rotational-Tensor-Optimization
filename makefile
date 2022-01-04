# --------------------------------------------------------------------------------------------------
# Makefile to build the tensor optimization program
#
# Author:    Johannes Gebert <gebert@hlrs.de>
# Date:      06.03.2021
# Last edit: 06.03.2021
#
# For use of make visit: https://www.gnu.org/software/make/
# --------------------------------------------------------------------------------------------------
# Set Architecture to compile for
# julius  - A Whiskey Lake Notebook, 4 cores, 16Gb memory, APU
trgt_arch = "julius"
# --------------------------------------------------------------------------------------------------
# Directories
mod_dir   = $(CURDIR)/mod/
obj_dir   = $(CURDIR)/obj/
lib_dir   = $(CURDIR)/lib/
bin_dir   = $(CURDIR)/bin/
f-src_dir = $(CURDIR)/f-src/
ext_f-src = $(CURDIR)/ext_f-src/
# --------------------------------------------------------------------------------------------------
# File extensions and suffixes
mod_ext   = .mod
obj_ext   = .o
sho_ext   = .so
f90_ext   = .f90
bin_suf   = _x86_64
# --------------------------------------------------------------------------------------------------
clean_cmd = rm -f
# --------------------------------------------------------------------------------------------------
# Compilers
#ifeq($(strip $(trgt_arch)) ,"julius" )
  compiler = "gfortran"
#endif
export compiler
# --------------------------------------------------------------------------------------------------
# Programming Environment - gnu, LLVM
PE = gnu
# --------------------------------------------------------------------------------------------------
# Compile flags GNU Compiler
ifeq ($(PE),gnu)
   c_flags_f90 = -J$(mod_dir) -I$(mod_dir) \
	               -g											   \
	               -O3										   \
	               -fbacktrace               \
                 -fbounds-check            \
                 -Wall #                   \				# Diagnoses
#	               -fdefault-integer-8       \				# incompatible with ISO_FORTRAN_ENV
#	               -fdefault-real-8          \
#	               -finstrument-functions    \
#	               -fopenmp
endif
# --------------------------------------------------------------------------------------------------
# Compile flags LLVM / Flang
# ifeq($(PE),LLVM)
#    c_Flags_f90 =
# endif
# --------------------------------------------------------------------------------------------------
# Executable
MAIN_bin = $(bin_dir)tensor_optimizer$(bin_suf)
# --------------------------------------------------------------------------------------------------
# Generate objects
#
f-objects = $(obj_dir)mod_standards$(obj_ext) \
			$(obj_dir)mod_stringmod$(obj_ext) \
			$(obj_dir)mod_file_routines$(obj_ext) \
			$(obj_dir)math_routines$(obj_ext) \
			$(obj_dir)opt_stiffness$(obj_ext) \
			$(obj_dir)tensor_optimizer$(obj_ext)
# --------------------------------------------------------------------------------------------------
# Begin Building
all: $(MAIN_bin)

# -------------------------------------------------------------------------------------------------
# Kinds Module
$(obj_dir)mod_standards$(obj_ext):$(f-src_dir)mod_standards$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $< "--------------------"
	$(compiler) $(c_flags_f90) -c $< -o $@
# --------------------------------------------------------------------------------------------------
# auxiliary routines module
$(obj_dir)mod_file_routines$(obj_ext):$(mod_dir)standards$(mod_ext) $(f-src_dir)mod_file_routines$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)mod_file_routines$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_file_routines$(f90_ext) -o $@
# --------------------------------------------------------------------------------------------------
# auxiliary routines module
$(obj_dir)math_routines$(obj_ext):$(mod_dir)standards$(mod_ext) $(f-src_dir)mod_math_routines$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)mod_math_routines$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_math_routines$(f90_ext) -o $@
#  -------------------------------------------------------------------------------------------------
# External source to parse input
$(obj_dir)mod_stringmod$(obj_ext):$(mod_dir)standards$(mod_ext)	$(f-src_dir)stringmod$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)stringmod$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)stringmod$(f90_ext) -o $@
# -------------------------------------------------------------------------------------------------
# inlining module
$(obj_dir)opt_stiffness$(obj_ext):$(mod_dir)standards$(mod_ext) \
								$(mod_dir)math_routines$(mod_ext)  \
								$(f-src_dir)mod_opt_stiffness$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)mod_opt_stiffness$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)mod_opt_stiffness$(f90_ext) -o $@

# --------------------------------------------------------------------------------------------------
# MAIN OBJECT
$(obj_dir)tensor_optimizer$(obj_ext):$(mod_dir)standards$(mod_ext) \
							$(mod_dir)math_routines$(mod_ext) \
							$(mod_dir)strings$(mod_ext) \
							$(mod_dir)file_routines$(mod_ext) \
							$(mod_dir)opt_stiffness$(mod_ext) \
							$(f-src_dir)tensor_optimizer$(f90_ext)
	@echo "---------------------------------------"
	@echo "-- Compiles: " $(f-src_dir)tensor_optimizer$(f90_ext) "--------------------"
	$(compiler) $(c_flags_f90) -c $(f-src_dir)tensor_optimizer$(f90_ext) -o $@
# ---------------------------------------------------------------------------------------------------
# Linking MAIN
$(MAIN_bin):$(f-objects)
	@echo "---------------------------------------"
	@echo "-- Linking MAIN binary"
	$(compiler) $(f-objects) -o $(MAIN_bin)

# ---------------------------------------------------------------------------------------------------
# Linking MAIN
	@echo "---------------------------------------"
	@echo "-- Successfully build all."
	@echo "---------------------------------------"


clean:
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning module directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(mod_dir)*$(mod_ext)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning object directory"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(f-objects)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning MAIN binary"
	@echo "----------------------------------------------------------------------------------"
	$(clean_cmd) $(MAIN_bin)
	@echo "----------------------------------------------------------------------------------"
	@echo "-- Cleaning completed."
	@echo "----------------------------------------------------------------------------------"
