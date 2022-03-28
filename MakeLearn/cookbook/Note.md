# From [makefile-cookbook](https://makefiletutorial.com/#makefile-cookbook)

## Get Started

### Syntax

```make
targets: prerequisites
    command
    command
    command
```

1. The targets are file names, separated by spaces.
2. The commands are a series of steps typically used to make the target(s).
3. The prerequisites are also file names, separated by spaces.

```make
blah: blah.o
    cc blah.o -o blah # Runs third

blah.o: blah.c
    cc -c blah.c -o blah.o # Runs second

blah.c:
    echo "int main() { return 0; }" > blah.c
```

### Variables

Variables can only be strings.   :=, = both works

```make
files := file1 file2
some_file: $(files)
    echo "Look at this variable: " $(files)
    touch some_file

file1:
    touch file1
file2:
    touch file2

clean:
    rm -f file1 file2 some_file
```

### Mutiple Targets

1. want multiple targets run

    ```make
    the: one two three

    one:
        touch one
    two:
        touch two
    three:
        touch three

    clean:
        rm -f one two three
    ```

2. multiple targets run follow the same rule

    ```make
    the: one two three

    one two three:
        touch $@

    clean:
        rm -f one two three
    ```

    $@ automatic variable that contains the target name

### Automatic Variables and Wildcards

#### Automatic Variables

```make
thing := one two three 

all: $(thing)
    echo $@
    echo $? # prerequisites newer than the target
    echo $< or $^ # all prerequisites
$(thing):
    touch $@

clean: 
    rm -f $(thing)
```

#### *

always wrap * with wildcard function

```make
thing := $(wildcard *.o)
all: $(thing)
    ls -l $(thing)
```

#### %

!!!!!!!!!!!!!!!!!!!!!!!!

### Fancy Rules

#### Implicit Rules

for C/C++ program

1. n.o is made automatically from n.c with a command of the form $(CC) -c $(CPPFLAGS) $(CFLAGS)
2. Compiling a C++ program: n.o is made automatically from n.cc or n.cpp with a command of the form $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)
3. Linking a single object file: n is made automatically from n.o by running the command $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)
4. CC: Program for compiling C programs; default cc
5. CXX: Program for compiling C++ programs; default g++
6. CFLAGS: Extra flags to give to the C compiler
7. CXXFLAGS: Extra flags to give to the C++ compiler
8. CPPFLAGS: Extra flags to give to the C preprocessor
9. LDFLAGS: Extra flags to give to compilers when they are supposed to invoke the linker

```make
CC = gcc
CFLAGS = -g

blah: blah.o
blah.o: blah.c
blah.c:
    echo "int main() { return 0; }" > $@

clean: 
    rm -f blah.*
```

#### Static Pattern Rules

```make
targets...: target-pattern: prereq-patterns
```

```make
CC = gcc
CFLAGS = -g

object = foo.o bar.o all.o
all : $(objects)

$(objects): %.o: %.c # foo.o mathches foo.c etc
                        # ! % usage

%.c:            # pattern rules
    touch $@
    echo $@

all.c :
    echo "int main() { return 0; }" > all.c

clean:
    rm -f foo* bar* all.c all
```

#### Static Pattern Rules and Filter

```make
obj_files = foo.result bar.o lose.o
src_files = foo.raw bar.c lose.c

all: $(obj_files)

$(filter %.result,$(obj_files)): %.result: %.raw
    touch $@
    echo "target: $@ prereq: $<" 

$(filter %.o,$(obj_files)): %.o: %.c
    touch $@
    echo "target: $@ prereq: $<"

%.c %.raw:
    touch $@

clean:
    rm -f $(src_files) $(obj_files)
```

### Commands and execution

#### Add an @ before a command to stop it from being printed

```make
all: 
    @echo "This make line will not be printed"
```

#### Each command is run in a new shell (or at least the effect is as such)

```make
all: 
    cd ..
    # The cd above does not affect this line, because each command is effectively run in a new shell
    echo `pwd`

    # This cd command affects the next because they are on the same line
    cd ..;echo `pwd`

    # Same as above
    cd ..; \
    echo `pwd`
```

#### Use export for recursive make

use the special $(MAKE) instead of make

```make
new_contents = "hello:\n\ttouch inside_file"
all:
    mkdir -p subdir
    printf $(new_contents) | sed -e 's/^ //' > subdir/makefile
    cd subdir && $(MAKE)
    @$(MAKE) clean

clean:
    rm -rf subdir
```

## Advanced

### Variable part 2

#### the difference of = and *= :*

```make
# Recursive variable. This will print "later" below
one = one ${later_variable}
# Simply expanded variable. This will not print "later" below
two := two ${later_variable}

later_variable = later

all: 
    echo $(one)
    echo $(two)
```

(using :=) allows you to append to a variable, we can also append with +=

```make
foo = start
foo += more
foo := $(foo) end
all: 
    echo $(foo)
```

?= only sets variables if they have not yet been set

```make
# Recursive variable. This will print "later" below
one = hello
one ?= will not be set
two ?= will be set

all: 
    echo $(one)
    echo $(two)
```

### Conditional part of Makefiles
