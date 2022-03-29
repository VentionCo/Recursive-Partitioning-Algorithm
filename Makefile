CC = emcc
OUT = -o dist/output.js 
OPTS =-s ALLOW_MEMORY_GROWTH=1 -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']" -s ENVIRONMENT='web' -s SINGLE_FILE=1

# Project name
PROJECT = program

SRCS := $(wildcard src/*.cpp)

$(PROJECT): build
		node build.js

# Targets
build: buildrepo
		$(CC) $(OUT) $(SRCS) $(OBJS)

clean:
		rm $(PROJECT) dist -Rf

buildrepo:
		mkdir -p dist