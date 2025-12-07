CC = gcc

# Flags de Compilação e Otimização
# Removi o DLIKWID_PERFMON daqui pois ele já está em LIKWID_FLAGS abaixo
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -Wall

# Configurações do LIKWID
# VERIFIQUE SE O CAMINHO ESTÁ CERTO NA SUA MÁQUINA
LIKWID_HOME = /home/soft/likwid
LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID_HOME)/include
LIKWID_LIBS = -L$(LIKWID_HOME)/lib -llikwid

# Bibliotecas de Linkagem (Math + Likwid)
LFLAGS = -lm $(LIKWID_LIBS)

PROG = cgSolver
MODULES = utils pcgc sislin
OBJS = $(addsuffix .o,$(MODULES)) $(PROG).o
SRCS = $(addsuffix .c,$(MODULES)) $(PROG).c $(addsuffix .h,$(MODULES))

# Arquivos para distribuição
DISTFILES = *.c *.h Makefile LEIAME.md benchmark.sh
DISTDIR = trabalho2_HPC

.PHONY: clean purge dist all debug

all: $(PROG)

# Regra genérica para gerar objetos
# Aqui unimos CFLAGS e LIKWID_FLAGS
%.o: %.c
	$(CC) -c $(CFLAGS) $(LIKWID_FLAGS) $< -o $@

# Regra de linkagem do executável
$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

# Target de debug
debug: CFLAGS = -O0 -g -Wall -DDEBUG
debug: $(PROG)

clean:
	@echo "Limpando objetos e arquivos temporários..."
	@rm -rf core *~ *.bak *.o

purge: clean
	@echo "Removendo executável..."
	@rm -f $(PROG)

dist: purge
	@echo "Gerando arquivo de distribuição ($(DISTDIR).tgz) ..."
	@mkdir -p $(DISTDIR)
	@cp $(DISTFILES) $(DISTDIR)
	@tar -czvf $(DISTDIR).tgz $(DISTDIR)
	@rm -rf $(DISTDIR)