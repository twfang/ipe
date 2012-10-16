#SHELL=/bin/bash
OBJSN=$(notdir $(OBJS))
OBJPN=$(notdir $(OBJP))
FC=$(FC90)
FFLAGS90=$(FFLAGS) 
LFLAGS=$(LFLAGS)
HOME_dev=$(HOME_dev)
all:
	echo eldyn compiler $(FC) $(FCparallel)
	echo eldyn1 par $(PAR)
	echo eldyn FFLAGS $(FFLAGS)
	echo building serial1
	echo $(OBJSN)
	echo building serial1
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJSN)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
ifeq ($(PAR),sms)
	echo building sms
	echo eldyn2 par $(PAR)
	echo $(FCparallel)
	make -f  $(HOME_dev)/Makefile.sms OBJ="$(OBJPN)" FC=$(FCparallel) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
else
	echo building serial2
	make -f  $(HOME_dev)/Makefile.original OBJ="$(OBJPN)" FC=$(FC) FFLAGS90="$(FFLAGS)" LFLAGS="$(LFLAGS)" HOME_dev=$(HOME_dev)
endif
