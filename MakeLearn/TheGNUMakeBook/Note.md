# GNU Make

*spend less time debugging makefiles and more time running fast builds.*

## An Introduction to Makefiles

## Writing Makefiles

### What Makefiles Contain

1. explicit rules
2. implicit rules
3. variable definitions
4. deriectives :instruction for make to do something special while reading the makefile
5. comments

### Including Other Makefiles

```make
include filenames
```

### Overriding Part of Another Makefile

A contains B, using a match-anything pattern rule to say than to remake any target that cannot be made from the information in A, make should look in B

```make
foo:
    frabnicate > foo

%:force
    @ $(MAKE) -f MakefileB $@

force : ;
```

1. ‘make foo’, make will find the rule in A
2. make bar, make looks in B
3. pattern rule has a pattern %
4. prerequisite force : to guarantee that the recipe will be run even if the target file already exists
5. force : ; prevent make use impilict rule to build it

### How make Reads a Makefile

## Writing Rules
