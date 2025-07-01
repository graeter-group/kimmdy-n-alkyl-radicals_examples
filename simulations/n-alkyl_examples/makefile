.PHONY: clean purge
clean:
	-rm ./systems/*/*__s*
	-rm ./systems/*/submit.sh

clean_all: clean
	-rm -r ./systems/*/slurm

purge: clean_all
	-rm ./systems/*/*.*
