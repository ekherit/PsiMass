
# следующие 3 строки лучше не трогать
ifeq (,$(WORKDIR))
WORKDIR = $(shell pushd ../ 1>/dev/null && pwd )
endif

#Всё что делается делается вслух :) - желательно раскоментировать
VERBOSEMAKE=""
#Если ключик -g при сборки не нужен - если уверен в себе
NDEBUG=""

#Компилятор для этого пакета
#CC_LOCAL := gcc
#Линкер для этого пакета
#LD_LOCAL := gcc

#PACKAGES_LOCAL=""

# Дополнительные опции для компилятора C/C++
COPTS  =   `root-config --cflags` -I./ 

# Дополнительные опции для компилятора Fortran 
FOPTS  =  -fno-second-underscore 
# если достали подчерки в фортране

# Дополнительные опции для линкера 
#LDOPTS = -lg2c $(shell $(ROOTSYS)/bin/root-config --libs --glibs)  -lMinuit -lNew
LDOPTS = -lg2c -Xlinker -rpath -Xlinker 'root-config --libdir'  -lMinuit -lNew
 #-shared #если хочется использовать разделяемую библиотеку
#Если необходимо добавить CERNLIB, то лучше воспользоваться этим ключиком
CERNLIBRARY = "" 
#Список библиотек, если вам не нравится стандартный набор (расширяем по требованию)
# jetset74 mathlib graflib geant321 grafX11 packlib
CERNLIBS = mathlib graflib grafX11 packlib 

#Если определена переменная ONLYBINARY, то библиотека в пакете отсутствует
#ONLYBINARY=""
# где вываливать исполняемые файлы
#BINDIR=./
#Раскомментируйте, если хочется собрать динамическую библиотеку
#LIB_SHARED = ""

# Дополнительные либы (вставляются после либ)
#LIB_LOCAL= 
# `cernlib jetset74 mathlib packlib` `root-config --libs`       

# Определим, какие программы мы будем собирать
# Здесь следует указывать только используемые при работе пакеты, 
# чтобы не засорять bin - все программки, которые пишутся дли примера
# и тестирования следует складывать в директорию example
BINARIES = 

# укажем, из каких модулей этого пакета они состоят
# (эти модули не будут включены в библиотеку)
# и какие библиотеки надо подключить при сборке
#templateTest_MODULES = test
#templateTest_LIBS    = KcTemplate



# следующую строку лучше не трогать
include $(WORKDIR)/KcReleaseTools/rules.mk








