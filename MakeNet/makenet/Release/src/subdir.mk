################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/AcceptPlacement.cpp \
../src/AcceptSparseNeighbourhood.cpp \
../src/QuadTree.cpp \
../src/QuadTreeNode.cpp \
../src/QuadTreeNodeData.cpp \
../src/manet.cpp 

CPP_DEPS += \
./src/AcceptPlacement.d \
./src/AcceptSparseNeighbourhood.d \
./src/QuadTree.d \
./src/QuadTreeNode.d \
./src/QuadTreeNodeData.d \
./src/manet.d 

OBJS += \
./src/AcceptPlacement.o \
./src/AcceptSparseNeighbourhood.o \
./src/QuadTree.o \
./src/QuadTreeNode.o \
./src/QuadTreeNodeData.o \
./src/manet.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp src/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-src

clean-src:
	-$(RM) ./src/AcceptPlacement.d ./src/AcceptPlacement.o ./src/AcceptSparseNeighbourhood.d ./src/AcceptSparseNeighbourhood.o ./src/QuadTree.d ./src/QuadTree.o ./src/QuadTreeNode.d ./src/QuadTreeNode.o ./src/QuadTreeNodeData.d ./src/QuadTreeNodeData.o ./src/manet.d ./src/manet.o

.PHONY: clean-src

