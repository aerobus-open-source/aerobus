FROM ubuntu:latest

RUN apt update
RUN apt install  -y clang g++ libomp-dev cppcheck cpplint doxygen texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra
RUN apt install -y build-essential
RUN apt install -y cmake
RUN apt install -y git

RUN mkdir -p /usr/share/github
WORKDIR /usr/share/github
RUN git clone https://github.com/google/benchmark.git
WORKDIR /usr/share/github/benchmark/
RUN mkdir -p build
WORKDIR /usr/share/github/benchmark/build/
RUN cmake .. -DCMAKE_BUILD_TYPE=RELEASE -DBENCHMARK_DOWNLOAD_DEPENDENCIES=ON
RUN make install
