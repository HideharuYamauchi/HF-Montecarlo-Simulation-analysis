#!bin/bash

GSL_INCLUDE_PATH=$(gsl-config --cflags)
ROOT_INCLUDE_PATH=$(root-config --cflags)

GSL_LIBRARY_PATH=$(gsl-config --libs)
ROOT_LIBRARY_PATH=$(root-config --libs)

g++ $(gsl-config --cflags) $(root-config --cflags) simulate.cc -o execute.out $(gsl-config --libs) $(root-config --libs)
echo "compile complete..."
echo "now execute..."
#g++ GSL_INCLUDE_PATH ROOT_INCLUDE_PATH simulate.cc -o execute.out GSL_LIBRARY_PATH ROOT_LIBRARY_PATH

#if [ $# -ne 6 ]                                                                                                                                                
#then
#    echo "usage: ./execute.out [TMmode] [PATH/to/magnetic_moment_dist] [PATH/to/Geant4_simulation_result] [PATH/to/Geant4_simulation_enviroment]"
#    exit 1
#fi

./execute.out 210 ../include/BRECON_MOM_20200716_6.txt /home/hideharu/position/result/pos/run0007.dat /home/hideharu/position/result/pos/environment_run0007.dat

# argv[3] = "path/to/run0000.dat"

# -o option is for set the name of binary
# プリプロセッサとは,コンパイルを行う前にソースコードに対して行われる前処理. c/c++では,ディレクティブ(# ...)という命令を処理する.
# -c option is for make the object file


# 静的ライブラリとはlibxxx.aのようなファイル. オブジェクトファイルの集合体で, arコマンドを使用して以下のように作成する.
# ar -r libsample.a a.o b.o c.o
# 静的ライブラリもオブジェクトファイルの集合体で,実行ファイルを作成するときに使用することができるが,コンパイル時に静的ライブラリがあるディレクトリを-Lオプションで指定する必要がある.


# a.c b.c libsample.a から実行ファイルmainを作るときは以下のようにしてディレクトリと静的ライブラリを指定する
# gcc a.c b.c -L. -lsample -o main
# 静的ライブラリは-lxxxの形で指定する. lib, .aの部分は省略できる.


# コンパイラがデフォルトで参照するパスは-v オプション(gcc -v xxx.c)を使用し,コンパイル処理の各段階において実行されたコマンドを出力することによって確認するのが手っ取り早い.
# 毎回ヘッダファイルの指定に#include "..."を使うのは面倒なため,その時に-I オプション(gcc -I../ -v xxx.c)を使用する.
