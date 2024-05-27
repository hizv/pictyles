CHARMDIR=~/bin/charm
CHARMC=$(CHARMDIR)/bin/charmc $(OPTS)
OPTS=-g -module CommonLBs -DNDEBUG -memory charmdebug 
#-O3 
# -tracemode projections -lz 
all: server.out

.PHONY: clean server.out

server.decl.h: server.ci server.hpp
	$(CHARMC) server.ci

server.def.h: server.ci server.hpp

# pic.decl.h: server.ci pic.hpp
# 	$(CHARMC) pic.ci

# pic.def.h: server.ci pic.hpp

server.o: server.cpp server.decl.h server.def.h server.hpp pic.hpp
	$(CHARMC) -c $<

server.out: server.o
	$(CHARMC) -language charm++ -o $@ server.o

run-server: server.out
	LD_LIBRARY_PATH=/home/sb56/anaconda3/envs/charm/lib ./charmrun +p4 ./server.out ++server ++server-port 10000

clean:
	rm -f *.decl.h *.def.h *.out *.o charmrun

clog:
	rm -f *.log.gz *.sts *.projrc
