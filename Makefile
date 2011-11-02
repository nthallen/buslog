LDFLAGS=-L/usr/local/lib
CPPFLAGS=-I/usr/local/include -I/usr/include/libxml2
LDLIBS=-lcurllog -lxml2 -lnort -lcurl
CXXFLAGS=-Wall -g

.PHONY : clean

buslog : buslog.cc

clean :
	rm buslog

