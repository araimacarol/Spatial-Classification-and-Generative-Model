################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../utility/myrandom.cpp \
../utility/parser.cpp \
../utility/token.cpp \
../utility/tokenlist.cpp 

CPP_DEPS += \
./utility/myrandom.d \
./utility/parser.d \
./utility/token.d \
./utility/tokenlist.d 

OBJS += \
./utility/myrandom.o \
./utility/parser.o \
./utility/token.o \
./utility/tokenlist.o 


# Each subdirectory must supply rules for building sources it contributes
utility/%.o: ../utility/%.cpp utility/subdir.mk
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -DDEBUGGING -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$@" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


clean: clean-utility

clean-utility:
	-$(RM) ./utility/myrandom.d ./utility/myrandom.o ./utility/parser.d ./utility/parser.o ./utility/token.d ./utility/token.o ./utility/tokenlist.d ./utility/tokenlist.o

.PHONY: clean-utility

