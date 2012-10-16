#SHELL=/bin/bash
OBJ=$(notdir $(OBJS))
FC=$(FC90)
FFLAGS90=$(FFLAGS) 
LFLAGS=$(LFLAGS)
HOME_dev=$(HOME_dev)
all:
	echo plasma compiler $(FC) $(FCparallel)
	echo plasma par $(PAR)
	echo plasma FFLAGS $(FFLAGS)
ifeq ($(PAR),sms)
	echo building sms
	echo $(FCparallel)
	make -f  $(HOME_dev)/Makefile.sms OBJ="$(OBJ)" FC=$(FCparallel) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
else
	echo building serial
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJ)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
endif
