TESTS = input_test source_test network_test rate_test

check_PROGRAMS = input_test source_test network_test rate_test

input_test_SOURCES = input_test.c

input_test_LDADD = ../src/input.o ../src/network.o

source_test_SOURCES = source_test.c

source_test_LDADD = ../src/input.o ../src/network.o

network_test_SOURCES = network_test.c

network_test_LDADD = ../src/input.o ../src/network.o

rate_test_SOURCES = rate_test.c

rate_test_LDADD = ../src/rates.o

CLEANFILES = input.ini source.mdl network.chm *.abun

