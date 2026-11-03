# Image for the paper job of .github/workflows/c-cpp.yml (and local reproduction of CI).
# Pinned base image and tool versions, so that rebuilding it gives the same environment.
FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential g++ clang libomp-dev cmake git ca-certificates \
        python3 python3-pip cppcheck cpplint \
        doxygen graphviz \
        texlive-latex-base texlive-fonts-recommended texlive-fonts-extra texlive-latex-extra \
        libsollya-dev \
    && rm -rf /var/lib/apt/lists/*

# Google Benchmark, pinned
RUN git clone --depth 1 --branch v1.9.1 https://github.com/google/benchmark.git /usr/share/github/benchmark \
    && cmake -S /usr/share/github/benchmark -B /usr/share/github/benchmark/build \
        -DCMAKE_BUILD_TYPE=Release -DBENCHMARK_ENABLE_TESTING=OFF \
    && cmake --build /usr/share/github/benchmark/build --parallel \
    && cmake --install /usr/share/github/benchmark/build \
    && rm -rf /usr/share/github/benchmark
