.PHONY : clean

FFTW_INC=../local/include
FFTW_LIB=../local/lib

SUNDIALS_INC=../local/include
SUNDIALS_LIB=../local/lib

CPPFLAGS = -fPIC -g -Iinclude/ \
           -I${FFTW_INC} -L${FFTW_LIB} \
           -I${SUNDIALS_INC} -L${SUNDIALS_LIB} \
           -lm \
           -lfftw3_omp -lfftw3 \
           -lsundials_cvodes -lsundials_nvecserial -lsundials_nvecopenmp \
           -lblas -llapack \
           -fopenmp

LDFLAGS = -shared

SOURCES = $(shell echo src/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

TARGET = libfidimag.so

all: $(TARGET)

clean:
	rm -rf $(OBJECTS) $(TARGET) *.dSYM

$(TARGET) : $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
