#Point to required headers
include_directories("Rmath/")

#Get all the C files
FILE(GLOB CFiles *.c)

#This file is optional
list(REMOVE_ITEM CFiles ${PROJECT_SOURCE_DIR}/src/utils/rmath/randmtzig.c)

# use 'for' loop initial declaration)
#set_property(TARGET Rmath PROPERTY C_STANDARD 99)
set (CMAKE_C_FLAGS "--std=gnu99 ${CMAKE_C_FLAGS}")

#Compile static library
add_library(rmath ${CFiles})
