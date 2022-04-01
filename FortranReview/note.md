# Fortran For Scients and Eengineer

vscode setup

1. install extension mordern fortran and fortran intellgence
2. pip install package fortran_language_server
3. add fortls path  for  fortran intellgence
4. vscode-modern-fortran-formatter and pip install --upgrade fprettify
5. add fprettify to your Fortran › Formatting: Path Specify the full path of where the formatter is installed shortcut shift+alt+f
6. goot to go

```bash
path=/work/home/xinyu/soft/XinYuEnv/bin/fortls
```

[terminal shortcuts](https://blog.csdn.net/chenh297/article/details/80076437)

## Basic Elements of Fortran

CONSTANTS AND VARIABLES

intrinsic data types:

interger real complexs logical character

integer

16bit 32bit 64bit use kinds to specify

character

If a character string must include an apostrophe, then that apostrophe may be represented by two consecutive single quotes.

```fortran
'man''s best friend'
```

Alternatively, the character string containing a single quote can be surrounded by double quotes. For example, the string “Man’s best friend”

```fortran
"Man's best friend"
```

named constants are created using the Parameter arrtibute of a type declaration statement

```fortran
type,Parameter::name = value [, name2 = value...]
```

### INTRINSIC FUNCTIONS

1. generic functions can be used with more than one type of input data.
2. specific functions can be used with only one type of input data.

### LIST-DIRECTED INPUT AND OUTPUT STATEMENTS

```fortran
READ (*,*) input_list
```

1. input_list is the list of variables into which the values being read are placed
2. first arterisk : the source of the input
3. second arterisk : the format of the input data

The term list-directed input means that the types of the variables in the variable list determine the required format of the input data

Each READ statement in a program begins reading from a new line of input data.

## Program Design and Branching Structure

### LOGICAL CONSTANTS, VARIABLES, AND OPERATORS

```fortran
.TRUE.
.FALSE.

/=  !not equal to
```
