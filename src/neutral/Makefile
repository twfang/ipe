#SHELL=/bin/bash
OBJS=$(OBJS)
OBJS=$(notdir $(OBJS))
OBJP=$(OBJP)
OBJP=$(notdir $(OBJP))
FC=$(FC90)
FFLAGS90=$(FFLAGS) 
LFLAGS=$(LFLAGS)
HOME_dev=$(HOME_dev)
all:
	echo neutral compiler $(FC) $(FCparallel)
	echo neutral par $(PAR)
	echo neutral FFLAGS $(FFLAGS)
	echo building serial
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJS)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
ifeq ($(PAR),sms)
	echo building sms
	echo $(FCparallel)
	make -f  $(HOME_dev)/Makefile.sms OBJ="$(OBJP)" FC=$(FCparallel) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
else
	echo building serial
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJP)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
endif
