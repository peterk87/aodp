# Builds manpages from:
#  - pod portions of perl utilities found in PERL_DIR
#  - pod correspondents to non-perl utilities found in POD_DIR


# version number; automatically generated from the latest git tag
# overwrites stuff between "" in the line below magic
# see the Makefile in the project root
# magic line below (do not modify)
# iMMmwyWGTTGeJ6TYQg2myA
VERSION= "2.5.0.1"

PERL_DIR= ../pl
POD_DIR= pod

PERL_SRC= $(wildcard $(PERL_DIR)/*)
POD_SRC= $(wildcard $(POD_DIR)/*)

PERL_TARGET= $(addsuffix .1.gz,$(notdir $(PERL_SRC)))
POD_TARGET=  $(addsuffix .1.gz,$(notdir $(POD_SRC)))

CPP_DIR= ../c
HELP_TARGET= $(CPP_DIR)/aodp.help.cpp

TARGET= $(PERL_TARGET) $(POD_TARGET) $(HELP_TARGET)

POD2MAN_ARGS= --center="Automated Design Pipeline (aodp) Manual" --release="aodp-$(VERSION)"

all: $(TARGET)

$(POD_TARGET): %.1.gz: $(POD_DIR)/%
	@pod2man $(POD2MAN_ARGS) $< | gzip > $@

$(PERL_TARGET): %.1.gz: $(PERL_DIR)/%
	@pod2man $(POD2MAN_ARGS) $< | gzip > $@

$(HELP_TARGET): $(CPP_DIR)/%.help.cpp: $(POD_DIR)/%
	@( echo "#include \"Application.h\"" ; echo ; echo const string Application::_help = ) > $@
	@pod2text $< | ../../util/lsre '"' '\"' | ../../util/lsre "$$" '\n"' | ../../util/lsre "^" '"' >> $@
	@echo ";" >> $@

clean:
	rm -f $(TARGET)
