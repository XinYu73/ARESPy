!output
#define MPI_PRINT(chara,value) if(parallel%isroot)print*,chara,value

!output time
#define MPI_TIME(chara,flag,cmd) CALL start_time(chara,falg); CALL cmd ; CALL end_time(chara,flag); CALL write_time(chara,flag)
#define MPI_SUMTIME(chara,flag,cmd) CALL start_time(chara,flag); CALL cmd ; CALL end_time(chara,flag)

!deallocate and allocate
#define RE_ALLOCATE(xx,size) if(allocated(xx))then; deallocate(xx); allocate(xx##size); endif

!output file and line
#define xltP print *, 'LINE',__LINE__,'file',__FILE__

!VERSION information
#define VERSION_GIT write(6,*) '>>>>>>>        69c6ffe 2022-03-15 09:08:21 +0800        <<<<<<<'
