# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.13

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2018.3.4\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2018.3.4\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\Anand\CLionProjects\Ising_Model

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/Ising_Model.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Ising_Model.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Ising_Model.dir/flags.make

CMakeFiles/Ising_Model.dir/main.f90.obj: CMakeFiles/Ising_Model.dir/flags.make
CMakeFiles/Ising_Model.dir/main.f90.obj: ../main.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object CMakeFiles/Ising_Model.dir/main.f90.obj"
	C:\MinGW\bin\gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c C:\Users\Anand\CLionProjects\Ising_Model\main.f90 -o CMakeFiles\Ising_Model.dir\main.f90.obj

CMakeFiles/Ising_Model.dir/main.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/Ising_Model.dir/main.f90.i"
	C:\MinGW\bin\gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E C:\Users\Anand\CLionProjects\Ising_Model\main.f90 > CMakeFiles\Ising_Model.dir\main.f90.i

CMakeFiles/Ising_Model.dir/main.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/Ising_Model.dir/main.f90.s"
	C:\MinGW\bin\gfortran.exe $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S C:\Users\Anand\CLionProjects\Ising_Model\main.f90 -o CMakeFiles\Ising_Model.dir\main.f90.s

# Object files for target Ising_Model
Ising_Model_OBJECTS = \
"CMakeFiles/Ising_Model.dir/main.f90.obj"

# External object files for target Ising_Model
Ising_Model_EXTERNAL_OBJECTS =

Ising_Model.exe: CMakeFiles/Ising_Model.dir/main.f90.obj
Ising_Model.exe: CMakeFiles/Ising_Model.dir/build.make
Ising_Model.exe: CMakeFiles/Ising_Model.dir/objects1.rsp
Ising_Model.exe: CMakeFiles/Ising_Model.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable Ising_Model.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\Ising_Model.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Ising_Model.dir/build: Ising_Model.exe

.PHONY : CMakeFiles/Ising_Model.dir/build

CMakeFiles/Ising_Model.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\Ising_Model.dir\cmake_clean.cmake
.PHONY : CMakeFiles/Ising_Model.dir/clean

CMakeFiles/Ising_Model.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\Anand\CLionProjects\Ising_Model C:\Users\Anand\CLionProjects\Ising_Model C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug C:\Users\Anand\CLionProjects\Ising_Model\cmake-build-debug\CMakeFiles\Ising_Model.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Ising_Model.dir/depend
