# :q:##############################################################
#
# Subdirectories for building
#
################################################################

error:
	@echo "+-----------------------------------------------------------------+"
	@echo "|                                                                 |"
	@echo "|              SMELIb and SMEEXE                                  |"
	@echo "|                                                                 |"
	@echo "| Usage: make all              install and test SMELib            |"
	@echo "|        make clean            clean *.o and test executables     |"
	@echo "|                                                                 |"
	@echo "|  Make sure the system-specific makefile.def has been edited     |"
	@echo "|  to reflect your system configuration.                          |"
	@echo "+-----------------------------------------------------------------+"

all: sme done

sme:
	cd ./sme_lib/src; make;
	cd ./bin; make clean; make;

clean:
	cd ./sme_lib/src; make clean;
	cd ./bin; make clean;

wipe:
	cd ./sme_lib/src; make wipe;
	cd ./bin; make wipe;

done:
	@echo "  "
	@echo " +---------------------------------------------------------------+"
	@echo " |                                                               |"
	@echo " |                   SMEEXE is built or installed.               |"
	@echo " |                                                               |"
	@echo " +---------------------------------------------------------------+"
	@echo "  "
	@echo "  "
