thing := one two three 

all: $(thing)
	echo $@
	echo $? # prerequisites newer than the target
	echo $^ # all prerequisites
$(thing):
	touch $@

clean: 
	rm -f $(thing)