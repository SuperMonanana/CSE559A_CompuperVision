OPT			= -O2 -w

CC			= g++

EXE1 	  		= main


LIBS                    = 

INCLUDE			=

CFLAGS 			= $(OPT) $(INCLUDE)

CPPS			= arguments.cpp array.cpp binaryfilereader.cpp image.cpp binaryfilewriter.cpp \
			commandlinearguments.cpp eigenfaces.cpp eigfaces.cpp error.cpp face.cpp faces.cpp \
			filereader.cpp filewriter.cpp filterkernel.cpp functions.cpp gaussian.cpp \
			imageformat.cpp jacob.cpp pyramid.cpp stdafx.cpp targa.cpp textfilereader.cpp \
			user.cpp users.cpp vector.cpp

OBJS			= $(CPPS:.cpp=.o)

all: $(OBJS) $(EXE1)

.cpp.o:
	$(CC) $(CFLAGS) -c $(@:.o=.cpp)

.c.o:
	$(CC) $(CFLAGS) -c $(@:.o=.c)

depend:
	makedepend -Y $(ARCH) *.cpp

$(EXE1): $(OBJS) $(EXE1).o
	$(CC) $(CFLAGS) -o $(EXE1) $(EXE1).o $(OBJS) $(LIBS) -lm 



clean: 
	rm -f core *.o $(EXE1) # DO NOT DELETE
# DO NOT DELETE

arguments.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
arguments.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
arguments.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
arguments.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
arguments.o: arguments.h commandlinearguments.h main.h
array.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
array.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
array.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
array.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
array.o: commandlinearguments.h main.h
binaryfilereader.o: stdafx.h functions.h error.h filereader.h
binaryfilereader.o: binaryfilereader.h textfilereader.h filewriter.h
binaryfilereader.o: binaryfilewriter.h array.h vector.h image.h face.h user.h
binaryfilereader.o: users.h faces.h eigfaces.h filterkernel.h imageformat.h
binaryfilereader.o: targa.h pyramid.h gaussian.h arguments.h
binaryfilereader.o: commandlinearguments.h main.h
binaryfilewriter.o: stdafx.h functions.h error.h filereader.h
binaryfilewriter.o: binaryfilereader.h textfilereader.h filewriter.h
binaryfilewriter.o: binaryfilewriter.h array.h vector.h image.h face.h user.h
binaryfilewriter.o: users.h faces.h eigfaces.h filterkernel.h imageformat.h
binaryfilewriter.o: targa.h pyramid.h gaussian.h arguments.h
binaryfilewriter.o: commandlinearguments.h main.h
commandlinearguments.o: stdafx.h functions.h error.h filereader.h
commandlinearguments.o: binaryfilereader.h textfilereader.h filewriter.h
commandlinearguments.o: binaryfilewriter.h array.h vector.h image.h face.h
commandlinearguments.o: user.h users.h faces.h eigfaces.h filterkernel.h
commandlinearguments.o: imageformat.h targa.h pyramid.h gaussian.h
commandlinearguments.o: arguments.h commandlinearguments.h main.h
eigenfaces.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
eigenfaces.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
eigenfaces.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
eigenfaces.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
eigenfaces.o: arguments.h commandlinearguments.h main.h
eigfaces.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
eigfaces.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
eigfaces.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
eigfaces.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
eigfaces.o: commandlinearguments.h main.h
error.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
error.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
error.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
error.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
error.o: commandlinearguments.h main.h
face.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
face.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
face.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
face.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
face.o: commandlinearguments.h main.h
faces.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
faces.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
faces.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
faces.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
faces.o: commandlinearguments.h main.h jacob.h
filereader.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
filereader.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
filereader.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
filereader.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
filereader.o: arguments.h commandlinearguments.h main.h
filewriter.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
filewriter.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
filewriter.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
filewriter.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
filewriter.o: arguments.h commandlinearguments.h main.h
filterkernel.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
filterkernel.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
filterkernel.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
filterkernel.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
filterkernel.o: arguments.h commandlinearguments.h main.h
functions.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
functions.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
functions.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
functions.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
functions.o: arguments.h commandlinearguments.h main.h
gaussian.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
gaussian.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
gaussian.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
gaussian.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
gaussian.o: commandlinearguments.h main.h
image.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
image.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
image.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
image.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
image.o: commandlinearguments.h main.h
imageformat.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
imageformat.o: textfilereader.h filewriter.h binaryfilewriter.h array.h
imageformat.o: vector.h image.h face.h user.h users.h faces.h eigfaces.h
imageformat.o: filterkernel.h imageformat.h targa.h pyramid.h gaussian.h
imageformat.o: arguments.h commandlinearguments.h main.h
jacob.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
jacob.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
jacob.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
jacob.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
jacob.o: commandlinearguments.h main.h jacob.h
main.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
main.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
main.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
main.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
main.o: commandlinearguments.h main.h
pyramid.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
pyramid.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
pyramid.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
pyramid.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
pyramid.o: commandlinearguments.h main.h
stdafx.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
stdafx.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
stdafx.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
stdafx.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
stdafx.o: commandlinearguments.h main.h
targa.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
targa.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
targa.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
targa.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
targa.o: commandlinearguments.h main.h
textfilereader.o: stdafx.h functions.h error.h filereader.h
textfilereader.o: binaryfilereader.h textfilereader.h filewriter.h
textfilereader.o: binaryfilewriter.h array.h vector.h image.h face.h user.h
textfilereader.o: users.h faces.h eigfaces.h filterkernel.h imageformat.h
textfilereader.o: targa.h pyramid.h gaussian.h arguments.h
textfilereader.o: commandlinearguments.h main.h
user.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
user.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
user.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
user.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
user.o: commandlinearguments.h main.h
users.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
users.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
users.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
users.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
users.o: commandlinearguments.h main.h
vector.o: stdafx.h functions.h error.h filereader.h binaryfilereader.h
vector.o: textfilereader.h filewriter.h binaryfilewriter.h array.h vector.h
vector.o: image.h face.h user.h users.h faces.h eigfaces.h filterkernel.h
vector.o: imageformat.h targa.h pyramid.h gaussian.h arguments.h
vector.o: commandlinearguments.h main.h
