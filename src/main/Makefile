#SHELL=/bin/bash
OBJS= $(OBJS)
OBJ=$(notdir $(OBJS))
FC=$(FC90)
FFLAGS90=$(FFLAGS) 
LFLAGS=$(LFLAGS)
HOME_dev=$(HOME_dev)
all:
	echo main compiler $(FC) $(FCparallel)
	echo main par $(PAR)
	echo main FFLAGS $(FFLAGS)
ifeq ($(PAR),sms)
	echo building sms
	echo $(FCparallel)
	make -f  $(HOME_dev)/Makefile.sms OBJ="$(OBJ)" FC=$(FCparallel) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
else
	echo building serial
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJS)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
endif