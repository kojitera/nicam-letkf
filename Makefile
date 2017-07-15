include ./Mkinclude

all:
	cd common && $(MAKE)
	pwd
	cd nicam/common && $(MAKE)
	pwd
	cd nicam/letkf && $(MAKE)
	pwd
	cd nicam/tool_letkf && ${MAKE}
	pwd
	cd nicam/obsope && ${MAKE}
	pwd

clean:
	cd common && $(MAKE) clean
	pwd
	cd nicam/common && $(MAKE) clean
	pwd
	cd nicam/letkf && $(MAKE) clean
	pwd
	cd nicam/tool_letkf && $(MAKE) clean
	pwd
	cd nicam/obsope && $(MAKE) clean
	pwd
	cd bin && rm *
	pwd
