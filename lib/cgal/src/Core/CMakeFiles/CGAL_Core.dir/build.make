# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/santhosh/Dropbox/Thesis/code/lib/cgal

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/santhosh/Dropbox/Thesis/code/lib/cgal

# Include any dependencies generated for this target.
include src/Core/CMakeFiles/CGAL_Core.dir/depend.make

# Include the progress variables for this target.
include src/Core/CMakeFiles/CGAL_Core.dir/progress.make

# Include the compile flags for this target's objects.
include src/Core/CMakeFiles/CGAL_Core.dir/flags.make

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o: src/Core/CMakeFiles/CGAL_Core.dir/flags.make
src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o: src/Core/all_files.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/santhosh/Dropbox/Thesis/code/lib/cgal/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o"
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/CGAL_Core.dir/all_files.cpp.o -c /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core/all_files.cpp

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CGAL_Core.dir/all_files.cpp.i"
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core/all_files.cpp > CMakeFiles/CGAL_Core.dir/all_files.cpp.i

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CGAL_Core.dir/all_files.cpp.s"
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core/all_files.cpp -o CMakeFiles/CGAL_Core.dir/all_files.cpp.s

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.requires:
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.requires

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.provides: src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.requires
	$(MAKE) -f src/Core/CMakeFiles/CGAL_Core.dir/build.make src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.provides.build
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.provides

src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.provides.build: src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o

# Object files for target CGAL_Core
CGAL_Core_OBJECTS = \
"CMakeFiles/CGAL_Core.dir/all_files.cpp.o"

# External object files for target CGAL_Core
CGAL_Core_EXTERNAL_OBJECTS =

lib/libCGAL_Core.so.10.0.0: src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o
lib/libCGAL_Core.so.10.0.0: /usr/lib/x86_64-linux-gnu/libmpfr.so
lib/libCGAL_Core.so.10.0.0: /usr/lib/x86_64-linux-gnu/libgmp.so
lib/libCGAL_Core.so.10.0.0: lib/libCGAL.so
lib/libCGAL_Core.so.10.0.0: /usr/lib/libboost_thread-mt.so
lib/libCGAL_Core.so.10.0.0: /usr/lib/libboost_system-mt.so
lib/libCGAL_Core.so.10.0.0: src/Core/CMakeFiles/CGAL_Core.dir/build.make
lib/libCGAL_Core.so.10.0.0: src/Core/CMakeFiles/CGAL_Core.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library ../../lib/libCGAL_Core.so"
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CGAL_Core.dir/link.txt --verbose=$(VERBOSE)
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && $(CMAKE_COMMAND) -E cmake_symlink_library ../../lib/libCGAL_Core.so.10.0.0 ../../lib/libCGAL_Core.so.10 ../../lib/libCGAL_Core.so

lib/libCGAL_Core.so.10: lib/libCGAL_Core.so.10.0.0

lib/libCGAL_Core.so: lib/libCGAL_Core.so.10.0.0

# Rule to build all files generated by this target.
src/Core/CMakeFiles/CGAL_Core.dir/build: lib/libCGAL_Core.so
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/build

src/Core/CMakeFiles/CGAL_Core.dir/requires: src/Core/CMakeFiles/CGAL_Core.dir/all_files.cpp.o.requires
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/requires

src/Core/CMakeFiles/CGAL_Core.dir/clean:
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core && $(CMAKE_COMMAND) -P CMakeFiles/CGAL_Core.dir/cmake_clean.cmake
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/clean

src/Core/CMakeFiles/CGAL_Core.dir/depend:
	cd /home/santhosh/Dropbox/Thesis/code/lib/cgal && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/santhosh/Dropbox/Thesis/code/lib/cgal /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/CGALCore /home/santhosh/Dropbox/Thesis/code/lib/cgal /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core /home/santhosh/Dropbox/Thesis/code/lib/cgal/src/Core/CMakeFiles/CGAL_Core.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/Core/CMakeFiles/CGAL_Core.dir/depend

