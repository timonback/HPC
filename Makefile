NAMES=BackBehrendt
MAILS=3back 3behrend

INF=informatik.uni-hamburg.de
COLOR=\e[0;32m
NC=\e[0m
CLEAN=find -mindepth 2 -name Makefile -execdir make clean \;
TAR=tar --exclude=*.swp --exclude=*.pdf --transform='s|$(DIR:/=)|$(NAMES)|' -czvf $(NAMES).tar.gz $(DIR)

all: cleansubdirs package send

cleansubdirs:
	@$(CLEAN) > /dev/null

help:
	@echo "Usage: make DIR=«Verzeichnis»"

package:
ifndef DIR
	$(error NUM is not defined. See make help)
endif

	@echo -e "$(COLOR)Creating Tar archive$(NC)"
	@$(TAR) > /dev/null 
	 
send:
	@echo -e "$(COLOR)Sending mail via thunderbird$(NC)"
	@thunderbird -compose "to=hr-abgabe@wr.$(INF),subject=$(NAMES),cc='$(foreach M,$(MAILS),$(M)@$(INF),)',attachment=$(shell pwd)/$(NAMES).tar.gz,body=Moin%2C%0Aanbei%20unsere%20Abgabe%0A%0AGruß%0A$(NAMES)"

clean:
	@echo -e "$(COLOR)Cleaning directory$(NC)"
	@$(RM) $(NAMES).tar.gz

