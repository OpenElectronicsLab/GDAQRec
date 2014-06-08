######################################################################
# Automatically generated by qmake (2.01a) Wed Apr 18 09:59:45 2007
######################################################################

TEMPLATE = app
TARGET = 
DEPENDPATH += .
INCLUDEPATH += .

# Input
HEADERS += plotter.h DAQReader.h
SOURCES += main.cpp plotter.cpp DAQReader.cpp
RESOURCES += plotter.qrc

# Input
HEADERS += DAQSettingsDialog/DAQSettingsDialog.h
FORMS += DAQSettingsDialog/DAQSettingsDialog.ui
SOURCES += DAQSettingsDialog/DAQSettingsDialog.cpp

isEmpty (DAQLIB) {
	DAQLIB+=comedi
}

unix:!macx {
	contains(DAQLIB, nidaqmxbase) {
		INCLUDEPATH += /usr/local/natinst/nidaqmxbase/include/ 
		LIBS += -lnidaqmxbase 
		DEFINES += USE_NIDAQMXBASE
	}
	
	contains(DAQLIB, comedi) {
		LIBS += -lcomedi
		DEFINES += USE_COMEDI
	}

	contains(DAQLIB, commandline) {
		DEFINES += USE_COMMANDLINE_DAQ
	}
}

macx {
	contains(DAQLIB, nidaqmxbase) {
		INCLUDEPATH += "/Applications/National\ Instruments/NI-DAQmx\ Base/includes/"
		LIBS += -framework nidaqmxbase
		LIBS += -framework nidaqmxbaselv
		DEFINES += USE_NIDAQMXBASE
	}
}	
