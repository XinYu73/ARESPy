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
    echo $^ # all prerequisites
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

