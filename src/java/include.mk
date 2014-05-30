#
# This is included by each package's specific Makefile.
#

JAVAHOME	    = /usr/local/java/1.3.1_18
JAVAC           = $(JAVAHOME)/bin/javac
CLASSPATH       = $(BASE_DIR)
DEPRECATION     = -deprecation
JFLAGS          = -d $(CLASSPATH) $(DEPRECATION) -classpath $(CLASSPATH)
DEBUG_OPTION    = -O
CP				= /bin/cp
JAR				= $(JAVAHOME)/bin/jar

.SUFFIXES: .java .class 

#$(CLASS_DIR)/*.class :
#	@mkdir -p $(CLASS_DIR)
#	$(JAVAC) $(JFLAGS) $(DEBUG_OPTION) *.java

$(CLASS_DIR)/%.class: %.java
	@mkdir -p $(CLASS_DIR)
	$(JAVAC) $(JFLAGS) $(DEBUG_OPTION) $?
