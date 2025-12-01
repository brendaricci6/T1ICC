CC = gcc

# Flags de Compilação e Otimização (Exigidas no Enunciado)
# -O3: Otimização máxima
# -march=native: Otimiza para a arquitetura da máquina onde está compilando
# -mavx: Habilita instruções vetoriais AVX
# -fopt-info-vec: Gera relatório sobre quais loops foram vetorizados (stderr)
CFLAGS = -O3 -march=native -mavx -fopt-info-vec -Wall

# Configurações do LIKWID
# Ajuste o caminho base (LIKWID_HOME) se necessário (ex: /usr)
LIKWID_HOME = /usr/local
LIKWID_FLAGS = -DLIKWID_PERFMON -I$(LIKWID_HOME)/include
LIKWID_LIBS = -L$(LIKWID_HOME)/lib -llikwid

# Bibliotecas de Linkagem (Math + Likwid)
LFLAGS = -lm $(LIKWID_LIBS)

PROG = cgSolver
MODULES = utils pcgc sislin
OBJS = $(addsuffix .o,$(MODULES)) $(PROG).o
# SRCS para dist
SRCS = $(addsuffix .c,$(MODULES)) $(PROG).c $(addsuffix .h,$(MODULES))

# Arquivos para distribuição
DISTFILES = *.c *.h Makefile LEIAME.md benchmark.sh
DISTDIR = trabalho2_HPC

.PHONY: clean purge dist all debug

all: $(PROG)

# Regra genérica para gerar objetos
%.o: %.c
	$(CC) -c $(CFLAGS) $(LIKWID_FLAGS) $< -o $@

# Regra de linkagem do executável
$(PROG): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)

# Target de debug (desativa otimizações, ativa símbolos e debug do código)
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