# Fortran For Scients and Eengineer

vscode setup

1. install extension mordern fortran and fortran intellgence
2. pip install package fortran_language_server
3. add fortls path  for  fortran intellgence
4. vscode-modern-fortran-formatter and pip install --upgrade fprettify
5. add fprettify to your Fortran › Formatting: Path Specify the full path of where the formatter is installed shortcut shift+alt+f
6. ctrl+tab 在边界组打开的文件中切换
7. alt left right 在文件和编辑位置之间切换
8. goot to go

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

```fortran
.and. !logical and
.or. 
.eqv. !logical equivalence
.neqv. logical nonequivalenca
.not.
```

#### Logical Values in Input and Output Statements

list-directed READ statement

If the input value is .TRUE., or the first character of the input value is T, then the logical variable will be set to .TRUE.

If the input value is .FALSE., or the first character of the input value is F, then the logical variable will be set to .FALSE.

list-directed WRITE statement

#### The SELECT CASE Construct

```fortran
    select case (int(a))
    case (:-1)
        WRITE (*, *) "It's below freezing today!"
    case (0)
        WRITE (*, *) "It's exactly at the freezing point."
    case (1:20)
        WRITE (*, *) "It's cool today."
    case (21:)
        WRITE (*, *) "It's hot today."
    case default
        write (*, *) "!^!"
    end select
```

1. case_value : case_value == case_expr
2. low_value : low_value<=case_expr
3. :high_value : case_expr<=high_value
4. a:b : a<= case_expr <= b

## Loops and Character Manipulation

```fortran
exit !break
cycle !continue
```

### character manipulation

concatenation //

1. ACHAR(ival) INT CHAR
2. IACHAR(char) CHAR INT
3. LEN(str1) CHAR INT
4. LEN_TRIM(str1) CHAR INT
5. TRIM(str1) CHAR CHAR return str1 with trailing blanks removed

## I/O basics

### Write   Barely read

```fortran
WRITE (*,100) i, result 100 
FORMAT (' The result for iteration ', I3,' is ', F7.3)


WRITE (*,100) i, x ! Format in FORMAT statement 
100 FORMAT (1X,I6,F10.2) 
CHARACTER(len=20) :: string ! Format in character variable 
string = '(1X,I6,F10.2)' 
WRITE (*,string) i, x 
WRITE (*,'(1X,I6,F10.2)') i, x ! Format in character constant
```

The format descriptor F7.3 specifies that a space seven characters wide should be used to print out the value of variable result, and that it should be printed with three digits to the right of the decimal point.

## Introduction to Arrays
