FROM ubuntu:24.04@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55 AS builder
ARG DEBIAN_FRONTEND=noninteractive
ARG LMCPP_REF=f6c6f51797912a369d6b27cf9e96f148a9a0e4d1
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates cmake g++ git ninja-build nlohmann-json3-dev python3 r-base-core=4.3.3-2build2 \
    && rm -rf /var/lib/apt/lists/*
RUN Rscript --version 2>&1 | grep -F '4.3.3'
RUN git clone https://github.com/lmctl-ai/lmcpp.git /src/lmcpp \
    && git -C /src/lmcpp checkout "$LMCPP_REF"
RUN cmake -S /src/lmcpp -B /build-lmcpp -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DLMCPP_BUILD_TESTS=OFF -DLMCPP_STATIC_STDLIB=ON -DCMAKE_INSTALL_PREFIX=/opt/lmcpp \
    && cmake --build /build-lmcpp --parallel 2 \
    && cmake --install /build-lmcpp
WORKDIR /src/rcli
COPY . .
RUN cmake -S . -B /build-rcli -G Ninja -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_PREFIX_PATH=/opt/lmcpp -DBUILD_TESTING=ON \
    && cmake --build /build-rcli --parallel 2

FROM ubuntu:24.04@sha256:534baea6a22c03a63003dbc8dbe78fe34bc0d7e595d9a9dc9834884ff530eb55 AS runtime
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends r-base-core=4.3.3-2build2 \
    && rm -rf /var/lib/apt/lists/* \
    && Rscript --version 2>&1 | grep -F '4.3.3'
COPY --from=builder /build-rcli/bin/rcli /usr/local/bin/rcli
ENTRYPOINT ["/usr/local/bin/rcli"]

FROM builder AS test
RUN ctest --test-dir /build-rcli --output-on-failure
CMD ["ctest", "--test-dir", "/build-rcli", "--output-on-failure"]
